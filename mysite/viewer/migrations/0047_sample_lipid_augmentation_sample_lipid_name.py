# Generated by Django 2.2.13 on 2020-08-20 01:43

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0046_remove_sample_lipid_sample_lipid_custom_augment'),
    ]

    operations = [
        migrations.AddField(
            model_name='sample_lipid_augmentation',
            name='sample_lipid_name',
            field=models.ForeignKey(default=0, on_delete=django.db.models.deletion.CASCADE, related_name='augment_sample_lipid_name', to='viewer.Sample_Lipid'),
            preserve_default=False,
        ),
    ]
