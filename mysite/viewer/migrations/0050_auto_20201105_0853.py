# Generated by Django 2.2.13 on 2020-11-05 13:53

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0049_auto_20200819_2207'),
    ]

    operations = [
        migrations.RenameField(
            model_name='sample_lipid_augmentation',
            old_name='tm_scattering_net_change',
            new_name='tmg_scattering_net_change',
        ),
    ]
