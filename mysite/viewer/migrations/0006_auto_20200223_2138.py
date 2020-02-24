# Generated by Django 2.2.7 on 2020-02-24 02:38

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0005_auto_20200223_2030'),
    ]

    operations = [
        migrations.AddField(
            model_name='symmetrical_parameters',
            name='separated',
            field=models.BooleanField(default=False, verbose_name='separated form factor'),
        ),
        migrations.AlterField(
            model_name='symmetrical_parameters',
            name='name',
            field=models.CharField(max_length=255, verbose_name='name'),
        ),
    ]