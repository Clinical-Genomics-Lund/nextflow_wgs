To run this workflow you need a valid Sentieon license.

To set up the required license server, follow the instructions as provided by Sentieon: https://support.sentieon.com/quick_start/

Then beyond this, the environment variable `SENTIEON_LICENSE` needs to be specified and set to the IP + port of the license server.

If you use the provided config files, this is specified as such within the config file:

```
env {
  SENTIEON_LICENSE = '10.139.0.101:8990'
}
```
