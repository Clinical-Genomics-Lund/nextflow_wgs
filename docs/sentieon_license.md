To run this workflow you need a valid Sentieon license. To set up the required license server, follow the instructions [as provided by Sentieon](https://support.sentieon.com/quick_start).

Then within your `nextflow.config`, the environment variable `SENTIEON_LICENSE` needs to be set to the IP + port of the license server.

```
env {
  SENTIEON_LICENSE = '10.139.0.101:8990'
}
```
