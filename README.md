# Figure of Merit Estimation Tool

Quickly estimate a figure of merit for given survey parameters.

## Usage

The program consists of three python scripts, which can either be run directly (requires only python) or using a webserver frontend (requires Docker). Usage for the two options is described below.

### Webserver frontend

The frontend is configured to run in a Docker container defined in the `webserver/Dockerfile` Dockerfile. Builld the image with
```
docker build -t baltay-fom-webserver .
```
and run with
```
docker run -dit --name baltay-fom -p 127.0.0.1:8080:8080 baltay-fom-webserver
```
The webserver can then be accessed at 127.0.0.1:8080 from your browser.

### Scripts

To run the scripts directly, without the frontend, the three python scripts
```
webserver/src/snflux.py
webserver/src/prism.py
webserver/src/snvar.py
```
must be executed in order. Parameters can be changed in the scripts, but are not accessible from the command line. Output is written to `snoutput.txt`.
