# Generated by Django 2.2.7 on 2019-12-27 18:54

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0010_auto_20191223_1524'),
    ]

    operations = [
        migrations.AlterField(
            model_name='data_lipid',
            name='data_set_title',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='data_lipid_set', to='viewer.Data_Set'),
        ),
    ]