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

        self.set_sample_lipids(sample_id)
        self.set_parameter(parameter_id)
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
        self.update_ranges()
        self.do_fit()

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
            self.xray_fig = plt.figure(figsize=(5.5,4.3))

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

            if parameter.fit_report_xray:
                self.fit_report_xray()
            
            xray_figures.append(mpld3.fig_to_html(xray_fig))
            plt.cla()

        # Neutron fit graphs
        for neutron_data in neutron_datas:
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
            if parameter.fit_report:
                self.fit_report_neutron()

            neutron_figures.append(mpld3.fig_to_html(neutron_fig))
            plt.cla()

        # Probability graphs
        prob_fig = plt.figure(figsize=(6,5))
        xray_sdp_graphs = []
        xray_sdp_data = {}
        neutron_sdp_graphs = []
        neutron_sdp_data = {}

        self.calculate_probability_graphs()
        prob_graph = mpld3.fig_to_html(prob_fig)


    def get_fit_main(self):
        # Fit data download
        if "fit_download" in request.POST:
            return self.fit_download(self.writer)
        elif "sdp_download" in request.POST:
            return self.sdp_download(self.writer)
        else:
            xray_graphs_and_forms = zip(xray_figures, xray_ranges, xray_scales, xray_datas)
            neutron_graphs_and_forms = zip(neutron_figures, neutron_ranges, neutron_scales, neutron_datas)

            plt.close('all')

            ## Done
            return render(request, 'viewer/fit_main.html', {
                'tutorial':tutorial,
                'xuser_tutorial':xuser_tutorial,
                'project':project,
                'sample':sample,
                'data_exists':data_exists,
                'parameter':parameter,
                'zero_parameter':zero_parameter,
                'parameter_update_form':parameter_update_form,
                'fit_result':fit_result,
                'show_stats':show_statistics,
                'show_probs':show_probabilities,
                'xray_graphs_and_forms':xray_graphs_and_forms,
                'neutron_graphs_and_forms':neutron_graphs_and_forms,
                'prob_graph':prob_graph,
                'xray_sdp_graphs':xray_sdp_graphs,
                'neutron_sdp_graphs':neutron_sdp_graphs,
                'additional_parameters':additional_parameters
            })

    def fit_download(self):
        # Filename
        file_name = str(project.project_title).replace(' ','-').replace(':','.')+'_FIT_download_'+'_'+str(sample.sample_title).replace(' ','-').replace(':','.')+'_'+str(parameter.name).replace(' ','-').replace(':','.')+now.strftime("%m-%d-%H.%M")+'.csv'

        # Create the HttpResponse object with the appropriate CSV header.
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename={0}'.format(file_name)

        writer = csv.writer(response)
        self.writer = writer
        writer.writerow(['VesicleViewer Fit output', now])
        writer.writerow(['Project Name', 'Sample Name', 'Parameter Set'])
        writer.writerow([project.project_title, sample.sample_title, parameter.name])
        writer.writerow([])

        if project.model_type == "SM":
            writer.writerow(['Calculated Parameters'])
            writer.writerow(['Db', additional_parameters[0]])
            writer.writerow(['2Dc', additional_parameters[1]])
            writer.writerow(['Dhh', additional_parameters[2]])
            writer.writerow(['Dh', parameter.headgroup_thickness])
            writer.writerow(['Al', parameter.lipid_area])
            writer.writerow([])
        elif project.model_type == "AS":
            writer.writerow(['Calculated Parameters'])
            writer.writerow([])

            writer.writerow(['Inner'])
            writer.writerow(['Db', additional_parameters[0]])
            writer.writerow(['2Dc', additional_parameters[1]])
            writer.writerow(['Dhh', additional_parameters[4]])
            writer.writerow(['Dh', parameter.in_headgroup_thickness])
            writer.writerow(['Al', parameter.in_lipid_area])
            writer.writerow([])

            writer.writerow(['Outter'])
            writer.writerow(['Db', additional_parameters[2]])
            writer.writerow(['2Dc', additional_parameters[3]])
            writer.writerow(['Dhh', additional_parameters[5]])
            writer.writerow(['Dh', parameter.out_headgroup_thickness])
            writer.writerow(['Al', parameter.out_lipid_area])
            writer.writerow([])

        writer.writerow({'Fit Statistics'})
        for line in parameter.fit_report:
            writer.writerow([line])

        writer.writerow([])

        writer.writerow(['Q', 'Experimental i', 'Experimental Error', 'Calculated i'])

        for xray_data in xray_datas:
            writer.writerow([xray_data.data_set_title])
            calculated_i_values = self.get_graph(parameter, sample_lipids, xray_data, project.system_tempurature, project.advanced_options)

            j = 0
            for i in range(xray_data.min_index, xray_data.max_index):
                writer.writerow([xray_data.q_value[i], xray_data.intensity_value[i], xray_data.error_value[i], calculated_i_values[j]])
                j = j+1

        for neutron_data in neutron_datas:
            writer.writerow([neutron_data.data_set_title])
            calculated_i_values = self.get_graph(parameter, sample_lipids, neutron_data, project.system_tempurature, project.advanced_options)

            j = 0
            for i in range(neutron_data.min_index, neutron_data.max_index):
                writer.writerow([neutron_data.q_value[i], neutron_data.intensity_value[i], neutron_data.error_value[i], calculated_i_values[j]])
                j = j+1

        return response

    def sdp_download(self):
        # Filename
        file_name = str(project.project_title).replace(' ','-').replace(':','.')+'_SDP_download_'+'_'+str(sample.sample_title).replace(' ','-').replace(':','.')+'_'+str(parameter.name).replace(' ','-').replace(':','.')+now.strftime("%m-%d-%H.%M")+'.csv'

        # Create the HttpResponse object with the appropriate CSV header.
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename={0}'.format(file_name)

        writer = csv.writer(response)
        self.writer = writer
        writer.writerow(['VesicleViewer SDP output', now])
        writer.writerow(['Project Name', 'Sample Name', 'Parameter Set'])
        writer.writerow([project.project_title, sample.sample_title, parameter.name])
        writer.writerow([])

        # Water probablilities
        writer.writerow(['Volume Probabilities'])
        self.water_probabilities(writer)

        # SDP and Scaled Probabilities
        writer.writerow([])
        writer.writerow(['Scattering Density Profile'])

        for xray_data in xray_datas:
            sdp_results = xray_sdp_data[xray_data]

            writer.writerow([])
            writer.writerow([xray_data.data_set_title])

            self.sdp_xray_data()

        for neutron_data in neutron_datas:
            sdp_results = neutron_sdp_data[neutron_data]

            writer.writerow([])
            writer.writerow([neutron_data.data_set_title])

            self.sdp_neutron_data()

        return response

    def update_ranges(self):
        xray_ranges = []
        xray_scales = []

        # Update q range for all x-ray datasets
        for xray_data in xray_datas:
            
            if xray_data.data_set_title in request.POST:
                xray_range_form = Data_Range_Form(request.POST)
                xray_scale_form = Data_Scale_Form(request.POST, instance=xray_data)
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

                    return redirect('viewer:fit_main', project_id=project.id, sample_id=sample.id, parameter_id=parameter.id)

            else:
                xray_range_form = Data_Range_Form()
                xray_scale_form = Data_Scale_Form(instance=xray_data)

            xray_ranges.append(xray_range_form)
            xray_scales.append(xray_scale_form)

        neutron_ranges = []
        neutron_scales = []

        # Update q range for all neutron datasets
        for neutron_data in neutron_datas:
            if neutron_data.data_set_title in request.POST:
                neutron_range_form = Data_Range_Form(request.POST)
                neutron_scale_form = Data_Scale_Form(request.POST, instance=neutron_data)
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

                    return redirect('viewer:fit_main', project_id=project.id, sample_id=sample.id, parameter_id=parameter.id)

            else:
                neutron_range_form = Data_Range_Form()
                neutron_scale_form = Data_Scale_Form(instance=neutron_data)

            neutron_ranges.append(neutron_range_form)
            neutron_scales.append(neutron_scale_form)

class SymmetricalFit(Fit):
    def set_sample_lipids(self):
        self.sample_lipids = Sample_Lipid.objects.filter(sample_title_id=sample_id, lipid_location='BOTH')
    
    def set_parameter(self, parameter_id):
        self.parameter = get_object_or_404(Symmetrical_Parameters, id=parameter_id)

    def set_zero_parameter(self):
        self.zero_parameter = True if \
            parameter.chain_volume == 0 or \
            parameter.headgroup_volume == 0 or \
            parameter.terminal_methyl_volume == 0 or \
            parameter.lipid_area == 0 or \
            parameter.sigma == 0 \
            else False

    def parameter_update(self):
        if "parameter_update" in self.request.POST:
            parameter_update_form = Symmetrical_Parameter_Fit_Form(self.request.POST, instance=self.parameter)
            if parameter_update_form.is_valid():
                self.parameter = parameter_update_form.save(commit=False)
                self.parameter.save()
        else:
            parameter_update_form = Symmetrical_Parameter_Fit_Form(instance=self.parameter)
    
    def do_fit(self):
        if "fit" in request.POST:
            # Do fit
            fit_result = symmetrical_fit(self.parameter, sample_lipids, datas, project.system_tempurature, project.advanced_options)
            fit_parameters = fit_result.params

            # Copy current instance
            new_parameter = deepcopy(self.parameter)

            # Set title
            new_parameter.name = now.strftime("%m/%d/%H:%M")

            # Set params
            new_parameter.terminal_methyl_volume = round(fit_parameters['terminal_methyl_volume'].value, 6)
            new_parameter.lipid_area = round(fit_parameters['area_per_lipid'].value, 6)
            new_parameter.headgroup_thickness = round(fit_parameters['headgroup_thickness'].value, 6)

            # Set report
            fit_report = lsq.fit_report(fit_result)
            new_parameter.fit_report = fit_report.split('\n')

            new_parameter.id = None
            new_parameter.save()

            for data in datas:
                data.scale = fit_parameters['scale_%i' % data.id].value
                data.background = fit_parameters['background_%i' % data.id].value

                data.save()

            # print(lsq.fit_report(fit_result))

            return redirect('viewer:fit_main', project_id=project.id, sample_id=sample.id, parameter_id=new_parameter.id)

    def fit_report_xray(self):
        plt.plot(
            xray_data.q_value[xray_data.min_index:xray_data.max_index],
            symmetrical_graph(self.parameter, sample_lipids, xray_data, project.system_tempurature, project.advanced_options),
            color='r',
            label='Best Fit',
            zorder=2
        )

    def fit_report_neutron(self):
        plt.plot(
            xray_data.q_value[xray_data.min_index:xray_data.max_index],
            symmetrical_graph(self.parameter, sample_lipids, xray_data, project.system_tempurature, project.advanced_options),
            color='r',
            label='Best Fit',
            zorder=2
        )

    def water_probabilities(self, writer):
        writer.writerow([])
        writer.writerow(['Headgroup'])
        writer.writerow(['z', 'Ph(z)'])
        for z, ph in zip (x_values, head_prob):
            writer.writerow([z, ph])

        writer.writerow([])
        writer.writerow(['Chains'])
        writer.writerow(['z', 'Phc(z)'])
        for z, pc in zip (x_values, chain_prob):
            writer.writerow([z, pc])

        writer.writerow([])
        writer.writerow(['Terminal Methyl'])
        writer.writerow(['z', 'Ptm(z)'])
        for z, ptm in zip (x_values, tm_prob):
            writer.writerow([z, ptm])

        writer.writerow([])
        writer.writerow(['Methylene'])
        writer.writerow(['z', 'Pch(z)'])
        for z, pm in zip (x_values, methylene_prob):
            writer.writerow([z, pm])

        writer.writerow([])
        writer.writerow(['Water'])
        writer.writerow(['z', 'Pw(z)'])
        for z, pw in zip (x_values, water_prob):
            writer.writerow([z, pw])

    def sdp_xray_data(self, writer):
        writer.writerow([])
        writer.writerow(['Combined SDP'])
        writer.writerow(['z', ''])
        for z, sdp in zip (x_values, sdp_results[0]):
            writer.writerow([z, sdp])

        writer.writerow([])
        writer.writerow(['Headgroup'])
        writer.writerow(['z', ''])
        for z, sdph in zip (x_values, sdp_results[1]):
            writer.writerow([z, sdph])

        writer.writerow([])
        writer.writerow(['Methylene'])
        writer.writerow(['z', ''])
        for z, sdpc in zip (x_values, sdp_results[2]):
            writer.writerow([z, sdpc])

        writer.writerow([])
        writer.writerow(['Terminal Methyl'])
        writer.writerow(['z', ''])
        for z, sdptm in zip (x_values, sdp_results[3]):
            writer.writerow([z, sdptm])

        writer.writerow([])
        writer.writerow(['Water'])
        writer.writerow(['z', ''])
        for z, sdpw in zip (x_values, sdp_results[4]):
            writer.writerow([z, sdpw])
    
    def sdp_neutron_data(self, writer):
        writer.writerow([])
        writer.writerow(['Combined SDP'])
        writer.writerow(['z', ''])
        for z, sdp in zip (x_values, sdp_results[0]):
            writer.writerow([z, sdp])

        writer.writerow([])
        writer.writerow(['Headgroup'])
        writer.writerow(['z', ''])
        for z, sdph in zip (x_values, sdp_results[1]):
            writer.writerow([z, sdph])

        writer.writerow([])
        writer.writerow(['Methylene'])
        writer.writerow(['z', ''])
        for z, sdpc in zip (x_values, sdp_results[2]):
            writer.writerow([z, sdpc])

        writer.writerow([])
        writer.writerow(['Terminal Methyl'])
        writer.writerow(['z', ''])
        for z, sdptm in zip (x_values, sdp_results[3]):
            writer.writerow([z, sdptm])

        writer.writerow([])
        writer.writerow(['Water'])
        writer.writerow(['z', ''])
        for z, sdpw in zip (x_values, sdp_results[4]):
            writer.writerow([z, sdpw])


class AsymmetricalFit(Fit):
    def set_sample_lipids(self):
        self.sample_lipids_in = Sample_Lipid.objects.filter(sample_title_id=sample_id, lipid_location='IN')
        self.sample_lipids_out = Sample_Lipid.objects.filter(sample_title_id=sample_id, lipid_location='OUT')
    
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
            parameter_update_form = Asymmetrical_Parameter_Fit_Form(self.request.POST, instance=self.parameter)
            if parameter_update_form.is_valid():
                self.parameter = parameter_update_form.save(commit=False)
                self.parameter.save()
        else:
            parameter_update_form = Asymmetrical_Parameter_Fit_Form(instance=parameter)

    def do_fit(self):
        if "fit" in request.POST:
            # Do fit
            fit_result = asymmetrical_fit(self.parameter, sample_lipids_in, sample_lipids_out, datas, project.system_tempurature, project.advanced_options)
            fit_parameters = fit_result.params

            # Copy current instance
            new_parameter = deepcopy(self.parameter)

            # Set title
            new_parameter.name = now.strftime("%m/%d/%H:%M")

            # Set params
            new_parameter.in_terminal_methyl_volume = round(fit_parameters['in_terminal_methyl_volume'].value, 6)
            new_parameter.in_lipid_area = round(fit_parameters['in_area_per_lipid'].value, 6)
            new_parameter.in_headgroup_thickness = round(fit_parameters['in_headgroup_thickness'].value, 6)

            new_parameter.out_terminal_methyl_volume = round(fit_parameters['out_terminal_methyl_volume'].value, 6)
            new_parameter.out_lipid_area = round(fit_parameters['out_area_per_lipid'].value, 6)
            new_parameter.out_headgroup_thickness = round(fit_parameters['out_headgroup_thickness'].value, 6)

            # Set report
            fit_report = lsq.fit_report(fit_result)
            new_parameter.fit_report = fit_report.split('\n')

            new_parameter.id = None
            new_parameter.save()

            for data in datas:
                data.scale = fit_parameters['scale_%i' % data.id].value
                data.background = fit_parameters['background_%i' % data.id].value

                data.save()

            # print(lsq.fit_report(fit_result))
            return redirect('viewer:fit_main', project_id=project.id, sample_id=sample.id, parameter_id=new_parameter.id)
    
    def fit_report_xray(self):
        plt.plot(
            xray_data.q_value[xray_data.min_index:xray_data.max_index],
            asymmetrical_graph(self.parameter, sample_lipids_in, sample_lipids_out, xray_data, project.system_tempurature, project.advanced_options),
            color='r',
            label='Best Fit',
            zorder=2
        )
    
    def fit_report_neutron(self):
        plt.plot(
            neutron_data.q_value[neutron_data.min_index:neutron_data.max_index],
            asymmetrical_graph(self.parameter, sample_lipids_in, sample_lipids_out, neutron_data, project.system_tempurature, project.advanced_options),
            color='r',
            label='Best Fit',
            zorder=2
        )

    def plot_fit_line(self):
        plt.plot(
            xray_data.q_value[xray_data.min_index:xray_data.max_index],
            asymmetrical_graph(self.parameter, sample_lipids_in, sample_lipids_out, xray_data, project.system_tempurature, project.advanced_options),
            color='r',
            label='Best Fit',
            zorder=2
        )

    def water_probabilities(self, writer):
        writer.writerow([])
        writer.writerow(['Headgroup'])
        writer.writerow(['z', 'Ph(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, ph in zip (in_x_values, in_head_prob):
            writer.writerow([z, ph])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, ph in zip (out_x_values, out_head_prob):
            writer.writerow([z, ph])

        writer.writerow([])
        writer.writerow(['Chains'])
        writer.writerow(['z', 'Phc(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, pc in zip (in_x_values, in_chain_prob):
            writer.writerow([z, pc])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, pc in zip (out_x_values, out_chain_prob):
            writer.writerow([z, pc])

        writer.writerow([])
        writer.writerow(['Terminal Methyl'])
        writer.writerow(['z', 'Ptm(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, ptm in zip (in_x_values, in_tm_prob):
            writer.writerow([z, ptm])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, ptm in zip (out_x_values, out_tm_prob):
            writer.writerow([z, ptm])

        writer.writerow([])
        writer.writerow(['Methylene'])
        writer.writerow(['z', 'Pch(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, pm in zip (in_x_values, in_methylene_prob):
            writer.writerow([z, pm])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, pm in zip (out_x_values, out_methylene_prob):
            writer.writerow([z, pm])

        writer.writerow([])
        writer.writerow(['Water'])
        writer.writerow(['z', 'Pw(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, pw in zip (in_x_values, in_water_prob):
            writer.writerow([z, pw])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, pw in zip (out_x_values, out_water_prob):
            writer.writerow([z, pw])
            
    def sdp_xray_data(self, writer):
        writer.writerow([])
        writer.writerow(['Combined SDP'])
        writer.writerow(['z', 'Pch(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdp in zip (in_x_values, sdp_results[0]):
            writer.writerow([z, sdp])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdp in zip (out_x_values, sdp_results[1]):
            writer.writerow([z, sdp])

        writer.writerow([])
        writer.writerow(['Headgroup'])
        writer.writerow(['z', 'Ph(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdph in zip (in_x_values, sdp_results[2]):
            writer.writerow([z, sdph])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdph in zip (out_x_values, sdp_results[3]):
            writer.writerow([z, sdph])

        writer.writerow([])
        writer.writerow(['Methylene'])
        writer.writerow(['z', 'Phc(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdpc in zip (in_x_values, sdp_results[4]):
            writer.writerow([z, sdpc])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdpc in zip (out_x_values, sdp_results[5]):
            writer.writerow([z, sdpc])

        writer.writerow([])
        writer.writerow(['Terminal Methyl'])
        writer.writerow(['z', 'Ptm(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdptm in zip (in_x_values, sdp_results[6]):
            writer.writerow([z, sdptm])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdptm in zip (out_x_values, sdp_results[7]):
            writer.writerow([z, sdptm])

        writer.writerow([])
        writer.writerow(['Water'])
        writer.writerow(['z', 'Pw(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdpw in zip (in_x_values, sdp_results[8]):
            writer.writerow([z, sdpw])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdpw in zip (out_x_values, sdp_results[9]):
            writer.writerow([z, sdpw])

    def sdp_neutron_data(self, writer):
        writer.writerow([])
        writer.writerow(['Combined SDP'])
        writer.writerow(['z', 'Pch(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdp in zip (in_x_values, sdp_results[0]):
            writer.writerow([z, sdp])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdp in zip (out_x_values, sdp_results[1]):
            writer.writerow([z, sdp])

        writer.writerow([])
        writer.writerow(['Headgroup'])
        writer.writerow(['z', 'Ph(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdph in zip (in_x_values, sdp_results[2]):
            writer.writerow([z, sdph])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdph in zip (out_x_values, sdp_results[3]):
            writer.writerow([z, sdph])

        writer.writerow([])
        writer.writerow(['Methylene'])
        writer.writerow(['z', 'Phc(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdpc in zip (in_x_values, sdp_results[4]):
            writer.writerow([z, sdpc])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdpc in zip (out_x_values, sdp_results[5]):
            writer.writerow([z, sdpc])

        writer.writerow([])
        writer.writerow(['Terminal Methyl'])
        writer.writerow(['z', 'Ptm(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdptm in zip (in_x_values, sdp_results[6]):
            writer.writerow([z, sdptm])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdptm in zip (out_x_values, sdp_results[7]):
            writer.writerow([z, sdptm])

        writer.writerow([])
        writer.writerow(['Water'])
        writer.writerow(['z', 'Pw(z)'])
        writer.writerow([])
        writer.writerow(['INNER'])
        for z, sdpw in zip (in_x_values, sdp_results[8]):
            writer.writerow([z, sdpw])
        writer.writerow([])
        writer.writerow(['OUTTER'])
        for z, sdpw in zip (out_x_values, sdp_results[9]):
            writer.writerow([z, sdpw])

# generator function
def generate_fit_main(request, project_id, sample_id, param_id):
    if project.model_type == "SM":
        fit = SymmetricalFit(request, project_id, sample_id, param_id)
    elif project.model_type == "AS":
        fit = AsymmetricalFit(request, project, sample, datas, param_id)
    
    return fit.get_fit_main()