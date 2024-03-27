import sys
import logging
import flask
import json
import pathlib
import snflux
import prism
import snvar

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

    try:
        snflux.snflux( flask.request.json )
    except Exception as ex:
        logger.exception( f"Failed to run snflux" )
        return { 'status': 'error', 'error': f'Failed to run snflux: {str(ex)}' }

    try:
        prism.prism( flask.request.json )
    except Exception as ex:
        logger.exception( f"Failed to run prism" )
        return { 'status': 'error', 'error': f'Failed to run prism: {str(ex)}' }

    sninvar = pathlib.Path( workdir / "sninvarim.txt" )
    sninvar.replace( workdir / "sninvar.txt" )

    try:
        snvar.snvar( flask.request.json )
    except Exception as ex:
        logger.exception( f"Failed to run snvar" )
        return { 'status': 'error', 'error': f'Failed to run snvar: {str(ex)}' }

    with open( workdir / "snoutput.txt", "r" ) as ifp:
        snoutput = ifp.read()

    return { 'status': 'ok', 'text': snoutput }
