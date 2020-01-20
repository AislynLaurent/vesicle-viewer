# Generated by Django 2.2.7 on 2020-01-14 23:31

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0029_auto_20200114_0032'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='lipid',
            name='total_volume',
        ),
        migrations.AddField(
            model_name='data_set',
            name='q_max_index',
            field=models.IntegerField(blank=True, null=True, verbose_name='q max'),
        ),
        migrations.AddField(
            model_name='data_set',
            name='q_min_index',
            field=models.IntegerField(blank=True, null=True, verbose_name='q min'),
        ),
    ]