# Generated by Django 3.0.3 on 2020-04-24 03:19

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('upload', '0012_auto_20200422_2103'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mechanism',
            name='pressure',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='mechanism',
            name='temperature',
            field=models.FloatField(blank=True, null=True),
        ),
    ]
