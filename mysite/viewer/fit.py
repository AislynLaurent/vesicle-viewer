from django.shortcuts import get_object_or_404, render, redirect
from django.http import HttpResponse
from django.utils import timezone

# note: might not need these?
from django.utils import text
from django.contrib import messages 

# Other libraries
import csv

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import mpld3
import numpy as np
import lmfit as lsq
from copy import deepcopy
import re

# note: test if you need these, importing * is annoying
# Models
from .models import *

# Forms
from .forms import *

# Other imports
from .symfit import *
from .asymfit import *
from .probabilities import *

class Fit:
    def __init__(self, request, project_id, sample_id, parameter_id):
        self.request = request

        ## Tutorials
        if self.request.user.is_anonymous:
            self.xuser_tutorial = False
        else:
            self.xuser = ExtendedUser.objects.get(user=self.request.user)
            self.xuser_tutorial = self.xuser.display_tutorial

        self.tutorial = True
        # Dismiss the tutorial
        if "dismiss_this" in self.request.POST:
            self.tutorial = False

        # Dismiss all tutorials
        if "dismiss_all" in self.request.POST:
            self.xuser.display_tutorial = False
            self.xuser.save()

        ## Overall
        # get initial project and sample
        self.project = get_object_or_404(Project, id=project_id)
        self.sample = get_object_or_404(Sample, id=sample_id)

        # Data
        self.datas = Data_Set.objects.filter(sample_title_id=sample_id)
        self.data_exists = True if self.datas else False

        self.set_parameter(parameter_id)
        self.set_sample_lipids()
        self.set_zero_parameter()

        self.data_exists = True if self.datas else False

        self.xray_datas = self.datas.filter(data_type='XR')
        self.neutron_datas = self.datas.filter(data_type='NU')

        # Declare
        self.now = timezone.now()
        self.fit_result = None
        self.show_statistics = False
        self.show_probabilities = False

        ### Forms
        self.parameter_update()

        # note: could improve this by running this in a separate function
        #       other than the constructor, and use the if statements for request.POST
        self.redirect = self.update_ranges()
        if self.redirect: return None

        self.redirect = self.do_fit()
        if self.redirect: return None

        # caluclated values
        self.calculated_i_values = []

        # Show stats / probabilities / graph
        if "statistics" in self.request.POST:
            self.show_statistics = True
            self.show_probabilities = False

        if "probabilities" in self.request.POST:
            self.show_probabilities = True
            self.show_statistics = False

        if "graphs" in self.request.POST:
            self.show_statistics = False
            self.show_probabilities = False

        ## GRAPHS
        self.xray_figures = []
        self.neutron_figures = []

        # X-Ray fit graphs
        for xray_data in self.xray_datas:
            xray_fig = plt.figure(figsize=(5.5,4.3))

            x = np.asarray(xray_data.q_value[xray_data.min_index:xray_data.max_index])
            y = np.asarray(xray_data.intensity_value[xray_data.min_index:xray_data.max_index])
            error = np.asarray(xray_data.error_value[xray_data.min_index:xray_data.max_index])

            # Data scatter plot
            plt.errorbar(
                x, 
                y,
                fmt='o',
                color='c',
                mfc='w',
                zorder=1
            )
            plt.errorbar(x, y-error, fmt='_', color='grey', zorder=0)
            plt.errorbar(x, y+error, fmt='_', color='grey', zorder=0)

            plt.xscale('log')
            plt.xlabel('q(A-1)')

            plt.yscale('log')
            plt.ylabel('Intensity (A.U.)')

            plt.title(xray_data.data_set_title)

            if self.parameter.fit_report:
                self.fit_report_xray(xray_data)
            
            self.xray_figures.append(mpld3.fig_to_html(xray_fig))
            plt.cla()

        # Neutron fit graphs
        for neutron_data in self.neutron_datas:
            neutron_fig = plt.figure(figsize=(5.5,4.3))

            x = np.asarray(neutron_data.q_value[neutron_data.min_index:neutron_data.max_index])
            y = np.asarray(neutron_data.intensity_value[neutron_data.min_index:neutron_data.max_index])
            error = np.asarray(neutron_data.error_value[neutron_data.min_index:neutron_data.max_index])

            # Data scatter plot
            plt.errorbar(
                x,
                y,
                fmt='o',
                color='c',
                mfc='w',
                zorder=1
            )
            plt.errorbar(x, y-error, fmt='_', color='grey', zorder=0)
            plt.errorbar(x, y+error, fmt='_', color='grey', zorder=0)

            plt.xscale('log')
            plt.xlabel('q(A-1)')

            plt.yscale('log')
            plt.ylabel('Intensity (A.U.)')

            plt.title(neutron_data.data_set_title)

            # Fit line
            if self.parameter.fit_report:
                self.fit_report_neutron()

            self.neutron_figures.append(mpld3.fig_to_html(neutron_fig))
            plt.cla()

        # Probability graphs
        prob_fig = plt.figure(figsize=(6,5))
        self.xray_sdp_graphs = []
        self.xray_sdp_data = {}
        self.neutron_sdp_graphs = []
        self.neutron_sdp_data = {}

        self.calculate_probability_graphs()
        self.prob_graph = mpld3.fig_to_html(prob_fig)


    def get_fit_main(self):
        # self.redirect could be made more modular (in its own function)
        if self.redirect: return redirect

        # Fit data download
        if "fit_download" in self.request.POST:
            return self.fit_download()
        elif "sdp_download" in self.request.POST:
            return self.sdp_download()
        else:
            self.xray_graphs_and_forms = zip(self.xray_figures, self.xray_ranges, self.xray_scales, self.xray_datas)
            self.neutron_graphs_and_forms = zip(self.neutron_figures, self.neutron_ranges, self.neutron_scales, self.neutron_datas)

            plt.close('all')

            ## Done
            return render(self.request, 'viewer/fit_main.html', {
                'tutorial':self.tutorial,
                'xuser_tutorial':self.xuser_tutorial,
                'project':self.project,
                'sample':self.sample,
                'data_exists':self.data_exists,
                'parameter':self.parameter,
                'zero_parameter':self.zero_parameter,
                'parameter_update_form':self.parameter_update_form,
                'fit_result':self.fit_result,
                'show_stats':self.show_statistics,
                'show_probs':self.show_probabilities,
                'xray_graphs_and_forms':self.xray_graphs_and_forms,
                'neutron_graphs_and_forms':self.neutron_graphs_and_forms,
                'prob_graph':self.prob_graph,
                'xray_sdp_graphs':self.xray_sdp_graphs,
                'neutron_sdp_graphs':self.neutron_sdp_graphs,
                'additional_parameters':self.additional_parameters
            })

    def fit_download(self):
        # Filename
        file_name = str(self.project.project_title).replace(' ','-').replace(':','.')+'_FIT_download_'+'_'+str(self.sample.sample_title).replace(' ','-').replace(':','.')+'_'+str(self.parameter.name).replace(' ','-').replace(':','.')+self.now.strftime("%m-%d-%H.%M")+'.csv'

        # Create the HttpResponse object with the appropriate CSV header.
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename={0}'.format(file_name)

        writer = csv.writer(response)
        writer.writerow(['VesicleViewer Fit output', self.now])
        writer.writerow(['Project Name', 'Sample Name', 'Parameter Set'])
        writer.writerow([self.project.project_title, self.sample.sample_title, self.parameter.name])
        writer.writerow([])

        if self.project.model_type == "SM":
            writer.writerow(['Calculated Parameters'])
            writer.writerow(['Db', self.additional_parameters[0]])
            writer.writerow(['2Dc', self.additional_parameters[1]])
            writer.writerow(['Dhh', self.additional_parameters[2]])
            writer.writerow(['Dh', self.parameter.headgroup_thickness])
            writer.writerow(['Al', self.parameter.lipid_area])
            writer.writerow([])
        elif self.project.model_type == "AS":
            writer.writerow(['Calculated Parameters'])
            writer.writerow([])

            writer.writerow(['Inner'])
            writer.writerow(['Db', self.additional_parameters[0]])
            writer.writerow(['2Dc', self.additional_parameters[1]])
            writer.writerow(['Dhh', self.additional_parameters[4]])
            writer.writerow(['Dh', self.parameter.in_headgroup_thickness])
            writer.writerow(['Al', self.parameter.in_lipid_area])
            writer.writerow([])

            writer.writerow(['Outter'])
            writer.writerow(['Db', self.additional_parameters[2]])
            writer.writerow(['2Dc', self.additional_parameters[3]])
            writer.writerow(['Dhh', self.additional_parameters[5]])
            writer.writerow(['Dh', self.parameter.out_headgroup_thickness])
            writer.writerow(['Al', self.parameter.out_lipid_area])
            writer.writerow([])

        writer.writerow({'Fit Statistics'})
        for line in self.parameter.fit_report:
            writer.writerow([line])

        writer.writerow([])

        writer.writerow(['Q', 'Experimental i', 'Experimental Error', 'Calculated i'])

        for xray_data in self.xray_datas:
            writer.writerow([xray_data.data_set_title])
            self.calculated_i_values = self.get_graph(self.parameter, self.sample_lipids, xray_data, self.project.system_tempurature, self.project.advanced_options)

            j = 0
            for i in range(xray_data.min_index, xray_data.max_index):
                writer.writerow([xray_data.q_value[i], xray_data.intensity_value[i], xray_data.error_value[i], self.calculated_i_values[j]])
                j = j+1

        for neutron_data in self.neutron_datas:
            writer.writerow([neutron_data.data_set_title])
            self.calculated_i_values = self.get_graph(self.parameter, self.sample_lipids, neutron_data, self.project.system_tempurature, self.project.advanced_options)

            j = 0
            for i in range(neutron_data.min_index, neutron_data.max_index):
                writer.writerow([neutron_data.q_value[i], neutron_data.intensity_value[i], neutron_data.error_value[i], self.calculated_i_values[j]])
                j = j+1

        return response

    def sdp_download(self):
        # Filename
        file_name = str(self.project.project_title).replace(' ','-').replace(':','.')+'_SDP_download_'+'_'+str(self.sample.sample_title).replace(' ','-').replace(':','.')+'_'+str(self.parameter.name).replace(' ','-').replace(':','.')+self.now.strftime("%m-%d-%H.%M")+'.csv'

        # Create the HttpResponse object with the appropriate CSV header.
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename={0}'.format(file_name)

        writer = csv.writer(response)
        writer.writerow(['VesicleViewer SDP output', self.now])
        writer.writerow(['Project Name', 'Sample Name', 'Parameter Set'])
        writer.writerow([self.project.project_title, self.sample.sample_title, self.parameter.name])
        writer.writerow([])

        # Water probablilities
        writer.writerow(['Volume Probabilities'])
        self.water_probabilities(writer)

        # SDP and Scaled Probabilities
        writer.writerow([])
        writer.writerow(['Scattering Density Profile'])

        for xray_data in self.xray_datas:
            self.sdp_results = self.xray_sdp_data[xray_data]

            writer.writerow([])
            writer.writerow([xray_data.data_set_title])

            self.sdp_xray_data()

        for neutron_data in self.neutron_datas:
            self.sdp_results = self.neutron_sdp_data[neutron_data]

            writer.writerow([])
            writer.writerow([neutron_data.data_set_title])

            self.sdp_neutron_data()

        return response

    def update_ranges(self):
        self.xray_ranges = []
        self.xray_scales = []

        # Update q range for all x-ray datasets
        for xray_data in self.xray_datas:
            
            if xray_data.data_set_title in self.request.POST:
                xray_range_form = Data_Range_Form(self.request.POST)
                xray_scale_form = Data_Scale_Form(self.request.POST, instance=xray_data)
                if xray_range_form.is_valid() and xray_scale_form.is_valid():
                    # Get scale values
                    scale_value = xray_scale_form.cleaned_data['scale']

                    if xray_data.scale_upperbound <= scale_value:
                        xray_data.scale_upperbound = (abs(scale_value)*1.5)
                    elif xray_data.scale_lowerbound >= scale_value:
                        xray_data.scale_lowerbound = -(abs(scale_value)*1.5)
                    
                    # Get range values
                    max_value = xray_range_form.cleaned_data['max_value']
                    min_value = xray_range_form.cleaned_data['min_value']

                    try:
                        # Find the indexes for the closest value in q_value
                        max_index = min(enumerate(xray_data.q_value), key=lambda x: abs(max_value - x[1]))
                        min_index = min(enumerate(xray_data.q_value), key=lambda x: abs(min_value - x[1]))

                        # Set the indexes in the db - +1 to take the closest on the inside for the low end
                        xray_data.max_index = max_index[0]
                        xray_data.min_index = min_index[0] + 1

                        xray_data.save()
                    except TypeError:
                        xray_data.save()

                    return redirect('viewer:fit_main', project_id=self.project.id, sample_id=self.sample.id, parameter_id=self.parameter.id)

            else:
                xray_range_form = Data_Range_Form()
                xray_scale_form = Data_Scale_Form(instance=xray_data)

            self.xray_ranges.append(xray_range_form)
            self.xray_scales.append(xray_scale_form)

        self.neutron_ranges = []
        self.neutron_scales = []

        # Update q range for all neutron datasets
        for neutron_data in self.neutron_datas:
            if neutron_data.data_set_title in self.request.POST:
                neutron_range_form = Data_Range_Form(self.request.POST)
                neutron_scale_form = Data_Scale_Form(self.request.POST, instance=neutron_data)
                if neutron_range_form.is_valid() and neutron_scale_form.is_valid():
                    # Get scale values
                    scale_value = neutron_scale_form.cleaned_data['scale']

                    if neutron_data.scale_upperbound <= scale_value:
                        neutron_data.scale_upperbound = (abs(scale_value)*1.5)
                    elif neutron_data.scale_lowerbound >= scale_value:
                        neutron_data.scale_lowerbound = -(abs(scale_value)*1.5)

                    # Get range values
                    max_value = neutron_range_form.cleaned_data['max_value']
                    min_value = neutron_range_form.cleaned_data['min_value']

                    try:
                        # Find the indexes for the closes value in q_value
                        max_index = min(enumerate(neutron_data.q_value), key=lambda x: abs(max_value - x[1]))
                        min_index = min(enumerate(neutron_data.q_value), key=lambda x: abs(min_value - x[1]))

                        # Set the indexes in the db - +1 to take the closest on the inside for the low end
                        neutron_data.max_index = max_index[0]
                        neutron_data.min_index = min_index[0] + 1

                        neutron_data.save()
                    except TypeError:
                        neutron_data.save()

                    return redirect('viewer:fit_main', project_id=self.project.id, sample_id=self.sample.id, parameter_id=self.parameter.id)

            else:
                neutron_range_form = Data_Range_Form()
                neutron_scale_form = Data_Scale_Form(instance=neutron_data)

            self.neutron_ranges.append(neutron_range_form)
            self.neutron_scales.append(neutron_scale_form)
        # no redirect if you reach here
        return None

class SymmetricalFit(Fit):
    def set_sample_lipids(self):
        self.sample_lipids = Sample_Lipid.objects.filter(sample_title_id=self.sample.id, lipid_location='BOTH')
    
    def set_parameter(self, parameter_id):
        self.parameter = get_object_or_404(Symmetrical_Parameters, id=parameter_id)

    def set_zero_parameter(self):
        self.zero_parameter = True if \
            self.parameter.chain_volume == 0 or \
            self.parameter.headgroup_volume == 0 or \
            self.parameter.terminal_methyl_volume == 0 or \
            self.parameter.lipid_area == 0 or \
            self.parameter.sigma == 0 \
            else False

    def parameter_update(self):
        if "parameter_update" in self.request.POST:
            self.parameter_update_form = Symmetrical_Parameter_Fit_Form(self.request.POST, instance=self.parameter)
            if self.parameter_update_form.is_valid():
                self.parameter = self.parameter_update_form.save(commit=False)
                self.parameter.save()
        else:
            self.parameter_update_form = Symmetrical_Parameter_Fit_Form(instance=self.parameter)
    
    def do_fit(self):
        if "fit" in self.request.POST:
            # Do fit
            self.fit_result = symmetrical_fit(self.parameter, self.sample_lipids, self.datas, self.project.system_tempurature, self.project.advanced_options)
            fit_parameters = self.fit_result.params

            # Copy current instance
            new_parameter = deepcopy(self.parameter)

            # Set title
            new_parameter.name = self.now.strftime("%m/%d/%H:%M")

            # Set params
            new_parameter.terminal_methyl_volume = round(fit_parameters['terminal_methyl_volume'].value, 6)
            new_parameter.lipid_area = round(fit_parameters['area_per_lipid'].value, 6)
            new_parameter.headgroup_thickness = round(fit_parameters['headgroup_thickness'].value, 6)

            # Set report
            fit_report = lsq.fit_report(self.fit_result)
            new_parameter.fit_report = fit_report.split('\n')

            new_parameter.id = None
            new_parameter.save()

            for data in self.datas:
                data.scale = fit_parameters['scale_%i' % data.id].value
                data.background = fit_parameters['background_%i' % data.id].value

                data.save()

            # print(lsq.fit_report(self.fit_result))

            return redirect('viewer:fit_main', project_id=self.project.id, sample_id=self.sample.id, parameter_id=new_parameter.id)

    def fit_report_xray(self, xray_data):
        plt.plot(
            xray_data.q_value[xray_data.min_index:xray_data.max_index],
            symmetrical_graph(self.parameter, self.sample_lipids, xray_data, self.project.system_tempurature, self.project.advanced_options),
            color='r',
            label='Best Fit',
            zorder=2
        )

    def fit_report_neutron(self, neutron_data):
        plt.plot(
            neutron_data.q_value[neutron_data.min_index:neutron_data.max_index],
            symmetrical_graph(parameter, sample_lipids, neutron_data, project.system_tempurature, project.advanced_options),
            color='r',
            label='Best Fit',
            zorder=2
        )
    
    def calculate_probability_graphs(self):
        self.x_values = np.arange(-40, 40, 0.2)
        
        # Calculate probabilities
        self.head_prob = head(
                self.parameter.chain_volume,
                self.parameter.headgroup_volume,
                self.parameter.lipid_area,
                self.parameter.headgroup_thickness,
                self.parameter.sigma,
                self.x_values
            )
        self.chain_prob = chain(
                self.parameter.chain_volume,
                self.parameter.lipid_area,
                self.parameter.sigma,
                self.x_values
            )
        self.tm_prob = terminal(
                self.parameter.terminal_methyl_volume,
                self.parameter.lipid_area,
                self.parameter.sigma,
                self.x_values
            )
        self.methylene_prob = methylene(
                self.parameter.chain_volume,
                self.parameter.terminal_methyl_volume,
                self.parameter.lipid_area,
                self.parameter.sigma,
                self.x_values
            )
        self.water_prob = water(
                self.parameter.chain_volume,
                self.parameter.headgroup_volume,
                self.parameter.lipid_area,
                self.parameter.headgroup_thickness,
                self.parameter.sigma,
                self.x_values
            )

        if self.zero_parameter:
            plt.plot()
            plt.title('! DIVIDE BY ZERO ERROR !')
        else:
            # headgroup
            plt.plot(
                self.x_values,
                self.head_prob,
                color='c',
                marker='.',
                markersize='5',
                label='Headgroup',
                zorder=0
            )
            # chain
            plt.plot(
                self.x_values,
                self.chain_prob,
                color='g',
                marker='v',
                markersize='5',
                label='Chains',
                zorder=1
            )
            # terminal methyl
            plt.plot(
                self.x_values,
                self.tm_prob,
                color='m',
                marker='s',
                markersize='5',
                label='Terminal Methyl',
                zorder=2
            )
            # methylene
            plt.plot(
                self.x_values,
                self.methylene_prob,
                color='k',
                marker='p',
                markersize='5',
                label='Methylene',
                zorder=3
            )
            # water
            plt.plot(
                self.x_values,
                self.water_prob,
                color='b',
                marker='x',
                markersize='5',
                label='Water',
                zorder=4
            )

            plt.legend(loc=1)
            plt.xlabel('Distance from bilayer center [Å]')
            plt.ylabel('Volume probability')

        for xray_data in self.xray_datas:
            xray_sdp_fig = plt.figure(figsize=(5.5,4.3))

            self.sdp_results = symmetrical_sdp(
                    self.parameter,
                    self.head_prob,
                    self.methylene_prob,
                    self.tm_prob,
                    self.water_prob,
                    self.sample_lipids,
                    xray_data,
                    self.project.system_tempurature,
                    self.project.advanced_options
                )

            additional_parameters = sym_additional_parameters(
                    self.parameter,
                    self.sample_lipids,
                    xray_data,
                    self.project.system_tempurature,
                    np.asarray(self.x_values),
                    np.asarray(self.head_prob),
                    self.project.advanced_options
                )

            self.xray_sdp_data[xray_data] = self.sdp_results

            if self.zero_parameter:
                plt.plot()
                plt.title(str(xray_data.data_set_title)+' ! DIVIDE BY ZERO ERROR !')
            else:
                plt.plot(
                    self.x_values,
                    self.sdp_results[0],
                    linewidth=4,
                    color='k',
                    markersize='5',
                    label='Combined SDP',
                    zorder=5
                )
                # headgroup
                plt.plot(
                    self.x_values,
                    self.sdp_results[1],
                    color='c',
                    marker='.',
                    markersize='5',
                    label='Headgroup',
                    zorder=0
                )
                # methyl
                plt.plot(
                    self.x_values,
                    self.sdp_results[2],
                    color='g',
                    marker='v',
                    markersize='5',
                    label='Methylene',
                    zorder=1
                )
                # terminal methyl
                plt.plot(
                    self.x_values,
                    self.sdp_results[3],
                    color='m',
                    marker='s',
                    markersize='5',
                    label='Terminal Methyl',
                    zorder=2
                )
                # water
                plt.plot(
                    self.x_values,
                    self.sdp_results[4],
                    color='b',
                    marker='x',
                    markersize='5',
                    label='Water',
                    zorder=4
                )

                plt.legend(loc=1)
                plt.xlabel('Distance from bilayer center [Å]')
                plt.ylabel('ED (e Å-3 )')

                plt.title(xray_data.data_set_title)

            self.xray_sdp_graphs.append(mpld3.fig_to_html(xray_sdp_fig))
            plt.cla()

        for neutron_data in neutron_datas:
            neutron_sdp_fig = plt.figure(figsize=(5.5,4.3))

            self.sdp_results = symmetrical_sdp(
                    self.parameter,
                    self.head_prob,
                    self.methylene_prob,
                    self.tm_prob,
                    self.water_prob,
                    sample_lipids,
                    neutron_data,
                    project.system_tempurature,
                    project.advanced_options
                )

            additional_parameters = sym_additional_parameters(
                    self.parameter,
                    sample_lipids,
                    neutron_data,
                    project.system_tempurature,
                    np.asarray(self.x_values),
                    np.asarray(self.head_prob),
                    project.advanced_options
                )

            neutron_sdp_data[neutron_data] = self.sdp_results

            if self.zero_parameter:
                plt.plot()
                plt.title(str(neutron_data.data_set_title)+' ! DIVIDE BY ZERO ERROR !')
            else:
                plt.plot(
                    self.x_values,
                    self.sdp_results[0],
                    linewidth=4,
                    color='k',
                    markersize='5',
                    label='Combined SDP',
                    zorder=5
                )
                # headgroup
                plt.plot(
                    self.x_values,
                    self.sdp_results[1],
                    color='c',
                    marker='.',
                    markersize='5',
                    label='Headgroup',
                    zorder=0
                )
                # methyl
                plt.plot(
                    self.x_values,
                    self.sdp_results[2],
                    color='g',
                    marker='v',
                    markersize='5',
                    label='Methylene',
                    zorder=1
                )
                # terminal methyl
                plt.plot(
                    self.x_values,
                    self.sdp_results[3],
                    color='m',
                    marker='s',
                    markersize='5',
                    label='Terminal Methyl',
                    zorder=2
                )
                # water
                plt.plot(
                    self.x_values,
                    self.sdp_results[4],
                    color='b',
                    marker='x',
                    markersize='5',
                    label='Water',
                    zorder=4
                )

                plt.legend(loc=1)
                plt.xlabel('Distance from bilayer center [Å]')
                plt.ylabel('NSLD (Å-3 x 10-5)')

                plt.title(neutron_data.data_set_title)

            neutron_sdp_graphs.append(mpld3.fig_to_html(neutron_sdp_fig))
            plt.cla()

    def water_probabilities(self, writer):
        writer.writerow([])
        writer.writerow(['Headgroup'])
        writer.writerow(['z', 'Ph(z)'])
        for z, ph in zip (self.x_values, self.head_prob):
            writer.writerow([z, ph])

        writer.writerow([])
        writer.writerow(['Chains'])
        writer.writerow(['z', 'Phc(z)'])
        for z, pc in zip (self.x_values, self.chain_prob):
            writer.writerow([z, pc])

        writer.writerow([])
        writer.writerow(['Terminal Methyl'])
        writer.writerow(['z', 'Ptm(z)'])
        for z, ptm in zip (self.x_values, self.tm_prob):
            writer.writerow([z, ptm])

        writer.writerow([])
        writer.writerow(['Methylene'])
        writer.writerow(['z', 'Pch(z)'])
        for z, pm in zip (self.x_values, self.methylene_prob):
            writer.writerow([z, pm])

        writer.writerow([])
        writer.writerow(['Water'])
        writer.writerow(['z', 'Pw(z)'])
        for z, pw in zip (self.x_values, self.water_prob):
            writer.writerow([z, pw])

    def sdp_xray_data(self, writer):
        writer.writerow([])
        writer.writerow(['Combined SDP'])
        writer.writerow(['z', ''])
        for z, sdp in zip (self.x_values, self.sdp_results[0]):
            writer.writerow([z, sdp])

        writer.writerow([])
        writer.writerow(['Headgroup'])
        writer.writerow(['z', ''])
        for z, sdph in zip (self.x_values, self.sdp_results[1]):
            writer.writerow([z, sdph])

        writer.writerow([])
        writer.writerow(['Methylene'])
        writer.writerow(['z', ''])
        for z, sdpc in zip (self.x_values, self.sdp_results[2]):
            writer.writerow([z, sdpc])

        writer.writerow([])
        writer.writerow(['Terminal Methyl'])
        writer.writerow(['z', ''])
        for z, sdptm in zip (self.x_values, self.sdp_results[3]):
            writer.writerow([z, sdptm])

        writer.writerow([])
        writer.writerow(['Water'])
        writer.writerow(['z', ''])
        for z, sdpw in zip (self.x_values, self.sdp_results[4]):
            writer.writerow([z, sdpw])
    
    def sdp_neutron_data(self, writer):
        writer.writerow([])
        writer.writerow(['Combined SDP'])
        writer.writerow(['z', ''])
        for z, sdp in zip (self.x_values, self.sdp_results[0]):
            writer.writerow([z, sdp])

        writer.writerow([])
        writer.writerow(['Headgroup'])
        writer.writerow(['z', ''])
        for z, sdph in zip (self.x_values, self.sdp_results[1]):
            writer.writerow([z, sdph])

        writer.writerow([])
        writer.writerow(['Methylene'])
        writer.writerow(['z', ''])
        for z, sdpc in zip (self.x_values, self.sdp_results[2]):
            writer.writerow([z, sdpc])

        writer.writerow([])
        writer.writerow(['Terminal Methyl'])
        writer.writerow(['z', ''])
        for z, sdptm in zip (self.x_values, self.sdp_results[3]):
            writer.writerow([z, sdptm])

        writer.writerow([])
        writer.writerow(['Water'])
        writer.writerow(['z', ''])
        for z, sdpw in zip (self.x_values, self.sdp_results[4]):
            writer.writerow([z, sdpw])


class AsymmetricalFit(Fit):
    def set_sample_lipids(self):
        self.sample_lipids_in = Sample_Lipid.objects.filter(sample_title_id=self.sample.id, lipid_location='IN')
        self.sample_lipids_out = Sample_Lipid.objects.filter(sample_title_id=self.sample.id, lipid_location='OUT')
    
    def set_parameter(self, parameter_id):
        self.parameter = get_object_or_404(Asymmetrical_Parameters, id=parameter_id)
    
    def set_zero_parameter(self):
        self.zero_parameter = True if \
            self.parameter.in_chain_volume == 0 or \
            self.parameter.in_headgroup_volume == 0 or \
            self.parameter.in_terminal_methyl_volume == 0 or \
            self.parameter.in_lipid_area == 0 or \
            self.parameter.out_chain_volume == 0  or \
            self.parameter.out_headgroup_volume == 0 or \
            self.parameter.out_terminal_methyl_volume == 0 or \
            self.parameter.out_lipid_area == 0 or \
            self.parameter.sigma == 0 \
            else False
    
    def parameter_update(self):
        if "parameter_update" in self.request.POST:
            self.parameter_update_form = Asymmetrical_Parameter_Fit_Form(self.request.POST, instance=self.parameter)
            if self.parameter_update_form.is_valid():
                self.parameter = self.parameter_update_form.save(commit=False)
                self.parameter.save()
        else:
            self.parameter_update_form = Asymmetrical_Parameter_Fit_Form(instance=self.parameter)

    def do_fit(self):
        if "fit" in self.request.POST:
            # Do fit
            self.fit_result = asymmetrical_fit(self.parameter, self.sample_lipids_in, self.sample_lipids_out, self.datas, self.project.system_tempurature, self.project.advanced_options)
            fit_parameters = self.fit_result.params

            # Copy current instance
            new_parameter = deepcopy(self.parameter)

            # Set title
            new_parameter.name = self.now.strftime("%m/%d/%H:%M")

            # Set params
            new_parameter.in_terminal_methyl_volume = round(fit_parameters['in_terminal_methyl_volume'].value, 6)
            new_parameter.in_lipid_area = round(fit_parameters['in_area_per_lipid'].value, 6)
            new_parameter.in_headgroup_thickness = round(fit_parameters['in_headgroup_thickness'].value, 6)

            new_parameter.out_terminal_methyl_volume = round(fit_parameters['out_terminal_methyl_volume'].value, 6)
            new_parameter.out_lipid_area = round(fit_parameters['out_area_per_lipid'].value, 6)
            new_parameter.out_headgroup_thickness = round(fit_parameters['out_headgroup_thickness'].value, 6)

            # Set report
            fit_report = lsq.fit_report(self.fit_result)
            new_parameter.fit_report = fit_report.split('\n')

            new_parameter.id = None
            new_parameter.save()

            for data in self.datas:
                data.scale = fit_parameters['scale_%i' % data.id].value
                data.background = fit_parameters['background_%i' % data.id].value

                data.save()

            # print(lsq.fit_report(self.fit_result))
            return redirect('viewer:fit_main', project_id=self.project.id, sample_id=self.sample.id, parameter_id=new_parameter.id)
    
    def fit_report_xray(self, xray_data):
        plt.plot(
            xray_data.q_value[xray_data.min_index:xray_data.max_index],
            asymmetrical_graph(self.parameter, self.sample_lipids_in, self.sample_lipids_out, xray_data, self.project.system_tempurature, self.project.advanced_options),
            color='r',
            label='Best Fit',
            zorder=2
        )
    
    def fit_report_neutron(self, neutron_data):
        plt.plot(
            neutron_data.q_value[neutron_data.min_index:neutron_data.max_index],
            asymmetrical_graph(self.parameter, self.sample_lipids_in, self.sample_lipids_out, neutron_data, self.project.system_tempurature, self.project.advanced_options),
            color='r',
            label='Best Fit',
            zorder=2
        )

    def plot_fit_line(self):
        plt.plot(
            xray_data.q_value[xray_data.min_index:xray_data.max_index],
            asymmetrical_graph(self.parameter, self.sample_lipids_in, self.sample_lipids_out, xray_data, self.project.system_tempurature, self.project.advanced_options),
            color='r',
            label='Best Fit',
            zorder=2
        )

    def calculate_probability_graphs(self):
        self.in_x_values = np.arange(-40, 0.2, 0.2)
        self.out_x_values = np.arange(-0.2, 40, 0.2)

        # Calculate probabilities
        self.in_head_prob = head(
                parameter.in_chain_volume,
                parameter.in_headgroup_volume,
                parameter.in_lipid_area,
                parameter.in_headgroup_thickness,
                parameter.sigma,
                self.in_x_values
            )
        self.out_head_prob = head(
                parameter.out_chain_volume,
                parameter.out_headgroup_volume,
                parameter.out_lipid_area,
                parameter.out_headgroup_thickness,
                parameter.sigma,
                self.out_x_values
            )
        self.in_chain_prob = chain(
                parameter.in_chain_volume,
                parameter.in_lipid_area,
                parameter.sigma,
                self.in_x_values
            )
        self.out_chain_prob = chain(
                parameter.out_chain_volume,
                parameter.out_lipid_area,
                parameter.sigma,
                self.out_x_values
            )
        self.in_tm_prob = terminal(
                parameter.in_terminal_methyl_volume,
                parameter.in_lipid_area,
                parameter.sigma,
                self.in_x_values
            )
        self.out_tm_prob = terminal(
                parameter.out_terminal_methyl_volume,
                parameter.out_lipid_area,
                parameter.sigma,
                self.out_x_values
            )
        self.in_methylene_prob = methylene(
                parameter.in_chain_volume,
                parameter.in_terminal_methyl_volume,
                parameter.in_lipid_area,
                parameter.sigma,
                self.in_x_values
            )
        self.out_methylene_prob = methylene(
                parameter.out_chain_volume,
                parameter.out_terminal_methyl_volume,
                parameter.out_lipid_area,
                parameter.sigma,
                self.out_x_values
            )
        self.in_water_prob = water(
                parameter.in_chain_volume,
                parameter.in_headgroup_volume,
                parameter.in_lipid_area,
                parameter.in_headgroup_thickness,
                parameter.sigma,
                self.in_x_values
            )
        self.out_water_prob = water(
                parameter.out_chain_volume,
                parameter.out_headgroup_volume,
                parameter.out_lipid_area,
                parameter.out_headgroup_thickness,
                parameter.sigma,
                self.out_x_values
            )

        if self.zero_parameter:
            plt.plot()
            plt.title('! DIVIDE BY ZERO ERROR !')
        else:
            # in headgroup
            plt.plot(
                self.in_x_values,
                self.in_head_prob,
                color='c',
                marker='.',
                markersize='5',
                label='Headgroup',
                zorder=0
            )
            # out headgroup
            plt.plot(
                self.out_x_values,
                self.out_head_prob,
                color='c',
                marker='.',
                markersize='5',
                zorder=0
            )
            # in chain
            plt.plot(
                self.in_x_values,
                self.in_chain_prob,
                color='g',
                marker='v',
                markersize='5',
                label='Chains',
                zorder=1
            )
            # out chain
            plt.plot(
                self.out_x_values,
                self.out_chain_prob,
                color='g',
                marker='v',
                markersize='5',
                zorder=1
            )
            # in terminal methyl
            plt.plot(
                self.in_x_values,
                self.in_tm_prob,
                color='m',
                marker='s',
                markersize='5',
                label='Terminal Methyl',
                zorder=2
            )
            # out terminal methyl
            plt.plot(
                self.out_x_values,
                self.out_tm_prob,
                color='m',
                marker='s',
                markersize='5',
                zorder=2
            )
            # in methylene
            plt.plot(
                self.in_x_values,
                self.in_methylene_prob,
                color='k',
                marker='p',
                markersize='5',
                label='Methylene',
                zorder=3
            )
            # out methylene
            plt.plot(
                self.out_x_values,
                self.out_methylene_prob,
                color='k',
                marker='p',
                markersize='5',
                zorder=3
            )
            # in water
            plt.plot(
                self.in_x_values,
                self.in_water_prob,
                color='b',
                marker='x',
                markersize='5',
                label='Water',
                zorder=4
            )
            # out water
            plt.plot(
                self.out_x_values,
                self.out_water_prob,
                color='b',
                marker='x',
                markersize='5',
                zorder=4
            )

            plt.legend(loc=1)
            plt.xlabel('Distance from bilayer center [Å]')
            plt.ylabel('Volume probability')

        for xray_data in xray_datas:
            xray_sdp_fig = plt.figure(figsize=(5.5,4.3))

            self.sdp_results = asymmetrical_sdp(
                    parameter,
                    self.in_head_prob,
                    self.in_methylene_prob,
                    self.in_tm_prob,
                    self.in_water_prob,
                    self.out_head_prob,
                    self.out_methylene_prob,
                    self.out_tm_prob,
                    self.out_water_prob,
                    self.sample_lipids_in,
                    self.sample_lipids_out, 
                    xray_data,
                    project.system_tempurature,
                    project.advanced_options
                )
            
            additional_parameters = asym_additional_parameters(
                    parameter,
                    self.sample_lipids_in,
                    self.sample_lipids_out, 
                    xray_data,
                    project.system_tempurature,
                    np.asarray(self.in_head_prob),
                    np.asarray(self.out_head_prob),
                    self.in_x_values,
                    self.out_x_values,
                    project.advanced_options
                )

            xray_sdp_data[xray_data] = self.sdp_results

            if zero_parameter:
                plt.plot()
                plt.title(str(xray_data.data_set_title)+' ! DIVIDE BY ZERO ERROR !')
            else:
                plt.plot(
                    self.in_x_values,
                    self.sdp_results[0],
                    linewidth=4,
                    color='k',
                    markersize='5',
                    label='Combined SDP',
                    zorder=5
                )
                plt.plot(
                    self.out_x_values,
                    self.sdp_results[1],
                    linewidth=4,
                    color='k',
                    markersize='5',
                    zorder=5
                )
                # headgroup
                plt.plot(
                    self.in_x_values,
                    self.sdp_results[2],
                    color='c',
                    marker='.',
                    markersize='5',
                    label='Headgroup',
                    zorder=0
                )
                plt.plot(
                    self.out_x_values,
                    self.sdp_results[3],
                    color='c',
                    marker='.',
                    markersize='5',
                    zorder=0
                )
                # methyl
                plt.plot(
                    self.in_x_values,
                    self.sdp_results[4],
                    color='g',
                    marker='v',
                    markersize='5',
                    label='Methylene',
                    zorder=1
                )
                plt.plot(
                    self.out_x_values,
                    self.sdp_results[5],
                    color='g',
                    marker='v',
                    markersize='5',
                    zorder=1
                )
                # terminal methyl
                plt.plot(
                    self.in_x_values,
                    self.sdp_results[6],
                    color='m',
                    marker='s',
                    markersize='5',
                    label='Terminal Methyl',
                    zorder=2
                )
                plt.plot(
                    self.out_x_values,
                    self.sdp_results[7],
                    color='m',
                    marker='s',
                    markersize='5',
                    zorder=2
                )
                # water
                plt.plot(
                    self.in_x_values,
                    self.sdp_results[8],
                    color='b',
                    marker='x',
                    markersize='5',
                    label='Water',
                    zorder=4
                )
                plt.plot(
                    self.out_x_values,
                    self.sdp_results[9],
                    color='b',
                    marker='x',
                    markersize='5',
                    zorder=4
                )

                plt.legend(loc=1)
                plt.xlabel('Distance from bilayer center [Å]')
                plt.ylabel('ED (e Å-3 )')

                plt.title(xray_data.data_set_title)

            self.xray_sdp_graphs.append(mpld3.fig_to_html(xray_sdp_fig))
            plt.cla()

        for neutron_data in neutron_datas:
            neutron_sdp_fig = plt.figure(figsize=(5.5,4.3))

            self.sdp_results = asymmetrical_sdp(
                    parameter,
                    self.in_head_prob,
                    self.in_methylene_prob,
                    self.in_tm_prob,
                    self.in_water_prob,
                    self.out_head_prob,
                    self.out_methylene_prob,
                    self.out_tm_prob,
                    self.out_water_prob,
                    self.sample_lipids_in,
                    self.sample_lipids_out, 
                    neutron_data,
                    project.system_tempurature,
                    project.advanced_options
                )

            additional_parameters = asym_additional_parameters(
                    parameter,
                    self.sample_lipids_in,
                    self.sample_lipids_out, 
                    neutron_data,
                    project.system_tempurature,
                    np.asarray(self.in_head_prob),
                    np.asarray(self.out_head_prob),
                    self.in_x_values,
                    self.out_x_values,
                    project.advanced_options
                )

            neutron_sdp_data[neutron_data] = self.sdp_results

            if zero_parameter:
                plt.plot()
                plt.title(str(neutron_data.data_set_title)+' ! DIVIDE BY ZERO ERROR !')
            else:
                plt.plot(
                    self.in_x_values,
                    self.sdp_results[0],
                    linewidth=4,
                    color='k',
                    markersize='5',
                    label='Combined SDP',
                    zorder=5
                )
                plt.plot(
                    self.out_x_values,
                    self.sdp_results[1],
                    linewidth=4,
                    color='k',
                    markersize='5',
                    zorder=5
                )
                # headgroup
                plt.plot(
                    self.in_x_values,
                    self.sdp_results[2],
                    color='c',
                    marker='.',
                    markersize='5',
                    label='Headgroup',
                    zorder=0
                )
                plt.plot(
                    self.out_x_values,
                    self.sdp_results[3],
                    color='c',
                    marker='.',
                    markersize='5',
                    zorder=0
                )
                # methyl
                plt.plot(
                    self.in_x_values,
                    self.sdp_results[4],
                    color='g',
                    marker='v',
                    markersize='5',
                    label='Methylene',
                    zorder=1
                )
                plt.plot(
                    self.out_x_values,
                    self.sdp_results[5],
                    color='g',
                    marker='v',
                    markersize='5',
                    zorder=1
                )
                # terminal methyl
                plt.plot(
                    self.in_x_values,
                    self.sdp_results[6],
                    color='m',
                    marker='s',
                    markersize='5',
                    label='Terminal Methyl',
                    zorder=2
                )
                plt.plot(
                    self.out_x_values,
                    self.sdp_results[7],
                    color='m',
                    marker='s',
                    markersize='5',
                    zorder=2
                )
                # water
                plt.plot(
                    self.in_x_values,
                    self.sdp_results[8],
                    color='b',
                    marker='x',
                    markersize='5',
                    label='Water',
                    zorder=4
                )
                plt.plot(
                    self.out_x_values,
                    self.sdp_results[9],
                    color='b',
                    marker='x',
                    markersize='5',
                    zorder=4
                )

                plt.legend(loc=1)
                plt.xlabel('Distance from bilayer center [Å]')
                plt.ylabel('NSLD (Å-3 x 10-5)')

                plt.title(neutron_data.data_set_title)

            neutron_sdp_graphs.append(mpld3.fig_to_html(neutron_sdp_fig))
            plt.cla()

    def water_probabilities(self, writer):
        writer.writerow([])
        writer.writerow(['Headgroup'])
        writer.writerow(['z', 'Ph(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, ph in zip (self.in_x_values, self.in_head_prob):
            writer.writerow([z, ph])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, ph in zip (self.out_x_values, self.out_head_prob):
            writer.writerow([z, ph])

        writer.writerow([])
        writer.writerow(['Chains'])
        writer.writerow(['z', 'Phc(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, pc in zip (self.in_x_values, self.in_chain_prob):
            writer.writerow([z, pc])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, pc in zip (self.out_x_values, self.out_chain_prob):
            writer.writerow([z, pc])

        writer.writerow([])
        writer.writerow(['Terminal Methyl'])
        writer.writerow(['z', 'Ptm(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, ptm in zip (self.in_x_values, self.in_tm_prob):
            writer.writerow([z, ptm])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, ptm in zip (self.out_x_values, self.out_tm_prob):
            writer.writerow([z, ptm])

        writer.writerow([])
        writer.writerow(['Methylene'])
        writer.writerow(['z', 'Pch(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, pm in zip (self.in_x_values, self.in_methylene_prob):
            writer.writerow([z, pm])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, pm in zip (self.out_x_values, self.out_methylene_prob):
            writer.writerow([z, pm])

        writer.writerow([])
        writer.writerow(['Water'])
        writer.writerow(['z', 'Pw(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, pw in zip (self.in_x_values, self.in_water_prob):
            writer.writerow([z, pw])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, pw in zip (self.out_x_values, self.out_water_prob):
            writer.writerow([z, pw])
            
    def sdp_xray_data(self, writer):
        writer.writerow([])
        writer.writerow(['Combined SDP'])
        writer.writerow(['z', 'Pch(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdp in zip (self.in_x_values, self.sdp_results[0]):
            writer.writerow([z, sdp])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdp in zip (self.out_x_values, self.sdp_results[1]):
            writer.writerow([z, sdp])

        writer.writerow([])
        writer.writerow(['Headgroup'])
        writer.writerow(['z', 'Ph(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdph in zip (self.in_x_values, self.sdp_results[2]):
            writer.writerow([z, sdph])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdph in zip (self.out_x_values, self.sdp_results[3]):
            writer.writerow([z, sdph])

        writer.writerow([])
        writer.writerow(['Methylene'])
        writer.writerow(['z', 'Phc(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdpc in zip (self.in_x_values, self.sdp_results[4]):
            writer.writerow([z, sdpc])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdpc in zip (self.out_x_values, self.sdp_results[5]):
            writer.writerow([z, sdpc])

        writer.writerow([])
        writer.writerow(['Terminal Methyl'])
        writer.writerow(['z', 'Ptm(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdptm in zip (self.in_x_values, self.sdp_results[6]):
            writer.writerow([z, sdptm])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdptm in zip (self.out_x_values, self.sdp_results[7]):
            writer.writerow([z, sdptm])

        writer.writerow([])
        writer.writerow(['Water'])
        writer.writerow(['z', 'Pw(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdpw in zip (self.in_x_values, self.sdp_results[8]):
            writer.writerow([z, sdpw])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdpw in zip (self.out_x_values, self.sdp_results[9]):
            writer.writerow([z, sdpw])

    def sdp_neutron_data(self, writer):
        writer.writerow([])
        writer.writerow(['Combined SDP'])
        writer.writerow(['z', 'Pch(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdp in zip (self.in_x_values, self.sdp_results[0]):
            writer.writerow([z, sdp])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdp in zip (self.out_x_values, self.sdp_results[1]):
            writer.writerow([z, sdp])

        writer.writerow([])
        writer.writerow(['Headgroup'])
        writer.writerow(['z', 'Ph(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdph in zip (self.in_x_values, self.sdp_results[2]):
            writer.writerow([z, sdph])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdph in zip (self.out_x_values, self.sdp_results[3]):
            writer.writerow([z, sdph])

        writer.writerow([])
        writer.writerow(['Methylene'])
        writer.writerow(['z', 'Phc(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdpc in zip (self.in_x_values, self.sdp_results[4]):
            writer.writerow([z, sdpc])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdpc in zip (self.out_x_values, self.sdp_results[5]):
            writer.writerow([z, sdpc])

        writer.writerow([])
        writer.writerow(['Terminal Methyl'])
        writer.writerow(['z', 'Ptm(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdptm in zip (self.in_x_values, self.sdp_results[6]):
            writer.writerow([z, sdptm])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdptm in zip (self.out_x_values, self.sdp_results[7]):
            writer.writerow([z, sdptm])

        writer.writerow([])
        writer.writerow(['Water'])
        writer.writerow(['z', 'Pw(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdpw in zip (self.in_x_values, self.sdp_results[8]):
            writer.writerow([z, sdpw])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdpw in zip (self.out_x_values, self.sdp_results[9]):
            writer.writerow([z, sdpw])

# generator function
def generate_fit_main(request, project_id, sample_id, param_id):
    if project.model_type == "SM":
        fit = SymmetricalFit(request, project_id, sample_id, param_id)
    elif project.model_type == "AS":
        fit = AsymmetricalFit(request, project, sample, datas, param_id)
    
    return fit.get_fit_main()