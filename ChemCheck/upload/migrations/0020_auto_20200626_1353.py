# Generated by Django 3.0.3 on 2020-06-26 13:53

from django.conf import settings
from django.db import migrations, models
import upload.models


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('upload', '0019_auto_20200616_1814'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mechanism',
            name='user',
            field=models.ForeignKey(default=1, editable=False, on_delete=models.SET(upload.models.get_sentinel_user), to=settings.AUTH_USER_MODEL),
        ),
    ]