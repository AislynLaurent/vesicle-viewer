# Generated by Django 2.2.7 on 2020-04-19 23:04

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0027_auto_20200415_1851'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='data_set',
            unique_together={('sample_title', 'data_set_title')},
        ),
    ]