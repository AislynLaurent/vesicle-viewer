# Generated by Django 2.2.7 on 2020-02-26 04:46

import django.contrib.postgres.fields
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0009_auto_20200224_1702'),
    ]

    operations = [
        migrations.AlterField(
            model_name='symmetrical_parameters',
            name='sample_title',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='sym_parameters', to='viewer.Sample'),
        ),
        migrations.CreateModel(
            name='Asymmetrical_Parameters',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255, verbose_name='name')),
                ('separated', models.BooleanField(default=False, verbose_name='separated form factor')),
                ('in_chain_volume', models.FloatField(default=0, verbose_name='in chain volume')),
                ('in_chain_volume_upperbound', models.FloatField(default=100, verbose_name='in cv upper bound')),
                ('in_chain_volume_lowerbound', models.FloatField(default=0, verbose_name='in cv lower bound')),
                ('in_chain_volume_lock', models.BooleanField(default=True, verbose_name='in cv lock')),
                ('in_headgroup_volume', models.FloatField(default=0, verbose_name='in headgroup volume')),
                ('in_headgroup_volume_upperbound', models.FloatField(default=100, verbose_name='in hv upper bound')),
                ('in_headgroup_volume_lowerbound', models.FloatField(default=0, verbose_name='in hv lower bound')),
                ('in_headgroup_volume_lock', models.BooleanField(default=True, verbose_name='in hv lock')),
                ('in_terminal_methyl_volume', models.FloatField(default=55, verbose_name='in terminal methyl volume')),
                ('in_terminal_methyl_volume_upperbound', models.FloatField(default=60, verbose_name='in tmv upper bound')),
                ('in_terminal_methyl_volume_lowerbound', models.FloatField(default=45, verbose_name='in tmv lower bound')),
                ('in_terminal_methyl_volume_lock', models.BooleanField(default=False, verbose_name='in tmv lock')),
                ('in_lipid_area', models.FloatField(default=60, verbose_name='in lipid area')),
                ('in_lipid_area_upperbound', models.FloatField(default=80, verbose_name='in la upper bound')),
                ('in_lipid_area_lowerbound', models.FloatField(default=40, verbose_name='in la lower bound')),
                ('in_lipid_area_lock', models.BooleanField(default=False, verbose_name='in la lock')),
                ('in_headgroup_thickness', models.FloatField(default=8, verbose_name='in headgroup thickness')),
                ('in_headgroup_thickness_upperbound', models.FloatField(default=15, verbose_name='in ht upper bound')),
                ('in_headgroup_thickness_lowerbound', models.FloatField(default=5, verbose_name='in ht lower bound')),
                ('in_headgroup_thickness_lock', models.BooleanField(default=False, verbose_name='in ht lock')),
                ('out_chain_volume', models.FloatField(default=0, verbose_name='out chain volume')),
                ('out_chain_volume_upperbound', models.FloatField(default=100, verbose_name='out cv upper bound')),
                ('out_chain_volume_lowerbound', models.FloatField(default=0, verbose_name='out cv lower bound')),
                ('out_chain_volume_lock', models.BooleanField(default=True, verbose_name='out cv lock')),
                ('out_headgroup_volume', models.FloatField(default=0, verbose_name='out headgroup volume')),
                ('out_headgroup_volume_upperbound', models.FloatField(default=100, verbose_name='out hv upper bound')),
                ('out_headgroup_volume_lowerbound', models.FloatField(default=0, verbose_name='out hv lower bound')),
                ('out_headgroup_volume_lock', models.BooleanField(default=True, verbose_name='out hv lock')),
                ('out_terminal_methyl_volume', models.FloatField(default=55, verbose_name='out terminal methyl volume')),
                ('out_terminal_methyl_volume_upperbound', models.FloatField(default=60, verbose_name='out tmv upper bound')),
                ('out_terminal_methyl_volume_lowerbound', models.FloatField(default=45, verbose_name='out tmv lower bound')),
                ('out_terminal_methyl_volume_lock', models.BooleanField(default=False, verbose_name='out tmv lock')),
                ('out_lipid_area', models.FloatField(default=60, verbose_name='out lipid area')),
                ('out_lipid_area_upperbound', models.FloatField(default=80, verbose_name='out la upper bound')),
                ('out_lipid_area_lowerbound', models.FloatField(default=40, verbose_name='out la lower bound')),
                ('out_lipid_area_lock', models.BooleanField(default=False, verbose_name='out la lock')),
                ('out_headgroup_thickness', models.FloatField(default=8, verbose_name='out headgroup thickness')),
                ('out_headgroup_thickness_upperbound', models.FloatField(default=15, verbose_name='out ht upper bound')),
                ('out_headgroup_thickness_lowerbound', models.FloatField(default=5, verbose_name='out ht lower bound')),
                ('out_headgroup_thickness_lock', models.BooleanField(default=False, verbose_name='out ht lock')),
                ('sigma', models.FloatField(default=2.5, verbose_name='sigma')),
                ('sigma_upperbound', models.FloatField(default=4, verbose_name='sig upper bound')),
                ('sigma_lowerbound', models.FloatField(default=2, verbose_name='sig lower bound')),
                ('sigma_lock', models.BooleanField(default=True, verbose_name='sig lock')),
                ('average_vesicle_radius', models.FloatField(default=500, verbose_name='average vesicle radius')),
                ('average_vesicle_radius_upperbound', models.FloatField(default=1000, verbose_name='avr upper bound')),
                ('average_vesicle_radius_lowerbound', models.FloatField(default=250, verbose_name='avr lower bound')),
                ('average_vesicle_radius_lock', models.BooleanField(default=True, verbose_name='avr lock')),
                ('relative_size', models.FloatField(default=0.25, verbose_name='relative size')),
                ('relative_size_upperbound', models.FloatField(default=0.7, verbose_name='rs upper bound')),
                ('relative_size_lowerbound', models.FloatField(default=0.1, verbose_name='rs lower bound')),
                ('relative_size_lock', models.BooleanField(default=True, verbose_name='rs lock')),
                ('fit_report', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(blank=True, max_length=500), blank=True, null=True, size=None, verbose_name='fitted parameter report')),
                ('sample_title', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='asym_parameters', to='viewer.Sample')),
            ],
            options={
                'verbose_name': 'parameter',
                'verbose_name_plural': 'parameters',
            },
        ),
    ]
