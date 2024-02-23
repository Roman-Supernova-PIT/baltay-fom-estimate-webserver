import sys
import logging
import flask
import json
import pathlib
import subprocess
import prism

workdir = pathlib.Path( __name__ ).resolve().parent

# "app" is the thing we give to gunicorn
app = flask.Flask( __name__, instance_relative_config=True )
app.logger.setLevel( logging.INFO )

@app.route( "/" )
def mainpage():
    return flask.render_template( 'fom.html' )


@app.route( "/calculate", methods=['POST'] )
def calculate( logger=logging.getLogger("main") ):
    response = {}

    failed = False
    try:
        result = subprocess.run( workdir / "snflux", capture_output=True )
    except Exception as ex:
        failed = True

    if failed or ( result.returncode != 0 ):
        logger.error( f"snflux returned {result.returncode}\n"
                      f"STDOUT\n------\n{result.stdout}\n"
                      f"STDERR\n------\n{result.stderr}" )
        response = { 'status': 'error', 'error': 'Error running snflux' };
        return response

    try:
        prism.prism( flask.request.json )
    except Exception as ex:
        logger.exception( f"Failed to run prism" )
        return { 'status': 'error', 'error': f'Failed to run prism: {str(ex)}' }

    sninvar = pathlib.Path( workdir / "sninvarim.txt" )
    sninvar.replace( workdir / "sninvar.txt" )

    failed = False
    try:
        result = subprocess.run( workdir / "snvar", capture_output=True )
    except Exception as ex:
        failed = True

    if failed or ( result.returncode != 0 ):
        logger.error( f"snvar returned {result.returncode}\n"
                      f"STDOUT\n------\n{result.stdout}\n"
                      f"STDERR\n------\n{result.stderr}\n" )
        response = { 'status': 'error', 'error': 'Error running snvar' }
        return response

    with open( workdir / "snoutput.txt", "r" ) as ifp:
        snoutput = ifp.read()

    return { 'status': 'ok', 'text': snoutput }
