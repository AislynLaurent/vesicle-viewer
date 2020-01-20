# Generated by Django 2.2.7 on 2020-01-19 22:39

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0031_auto_20200117_1224'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='data_set',
            name='q_ranged_values',
        ),
        migrations.RemoveField(
            model_name='molecule',
            name='total_volume',
        ),
        migrations.AddField(
            model_name='molecule',
            name='total_volume_equation',
            field=models.CharField(default='x', max_length=200, verbose_name='total volume equation'),
        ),
    ]