# Generated by Django 2.2.24 on 2022-01-30 02:11

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0051_auto_20211028_2217'),
    ]

    operations = [
        migrations.AddField(
            model_name='symmetrical_parameters',
            name='use_structure_factor',
            field=models.BooleanField(default=False, verbose_name='structure factor'),
        ),
    ]
