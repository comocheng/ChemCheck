# Generated by Django 3.0.3 on 2020-06-16 18:14

from django.conf import settings
from django.db import migrations, models
import upload.models


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('upload', '0018_mechanism_user'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mechanism',
            name='user',
            field=models.ForeignKey(null=True, on_delete=models.SET(upload.models.get_sentinel_user), to=settings.AUTH_USER_MODEL),
        ),
    ]
