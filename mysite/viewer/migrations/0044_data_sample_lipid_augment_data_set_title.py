# Generated by Django 2.2.13 on 2020-08-19 21:37

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0043_auto_20200819_1656'),
    ]

    operations = [
        migrations.AddField(
            model_name='data_sample_lipid_augment',
            name='data_set_title',
            field=models.ForeignKey(default=0, on_delete=django.db.models.deletion.CASCADE, related_name='data_set_augment_title', to='viewer.Data_Set'),
            preserve_default=False,
        ),
    ]
