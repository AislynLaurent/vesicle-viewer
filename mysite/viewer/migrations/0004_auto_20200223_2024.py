# Generated by Django 2.2.7 on 2020-02-24 01:24

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0003_auto_20200223_2023'),
    ]

    operations = [
        migrations.AlterField(
            model_name='symmetrical_parameters',
            name='average_vesicle_radius_lock',
            field=models.BooleanField(blank=True, default=True, null=True, verbose_name='avr lock'),
        ),
        migrations.AlterField(
            model_name='symmetrical_parameters',
            name='average_vesicle_radius_lowerbound',
            field=models.FloatField(blank=True, default=-1, null=True, verbose_name='avr lower bound'),
        ),
        migrations.AlterField(
            model_name='symmetrical_parameters',
            name='average_vesicle_radius_upperbound',
            field=models.FloatField(blank=True, default=1, null=True, verbose_name='avr upper bound'),
        ),
        migrations.AlterField(
            model_name='symmetrical_parameters',
            name='relative_size',
            field=models.FloatField(blank=True, default=0, null=True, verbose_name='relative size'),
        ),
        migrations.AlterField(
            model_name='symmetrical_parameters',
            name='relative_size_lock',
            field=models.BooleanField(blank=True, default=True, null=True, verbose_name='rs lock'),
        ),
        migrations.AlterField(
            model_name='symmetrical_parameters',
            name='relative_size_lowerbound',
            field=models.FloatField(blank=True, default=-1, null=True, verbose_name='rs lower bound'),
        ),
        migrations.AlterField(
            model_name='symmetrical_parameters',
            name='relative_size_upperbound',
            field=models.FloatField(blank=True, default=1, null=True, verbose_name='rs upper bound'),
        ),
    ]
