# Generated by Django 2.2.7 on 2019-12-18 18:09

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='symmetrical_parameters',
            name='headgroup_thickness',
            field=models.FloatField(default=6, verbose_name='headgroup thickness'),
        ),
        migrations.AlterField(
            model_name='symmetrical_parameters',
            name='lipid_area',
            field=models.FloatField(default=60, verbose_name='lipid area'),
        ),
        migrations.AlterField(
            model_name='symmetrical_parameters',
            name='scale',
            field=models.FloatField(default=1, verbose_name='scale'),
        ),
    ]