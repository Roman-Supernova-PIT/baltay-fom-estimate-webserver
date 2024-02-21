import sys
import logging
import flask
import json
import pathlib
import subprocess

workdir = pathlib.Path( __name__ ).resolve().parent
# sys.path.append( workdir )
import prism

app = flask.Flask( __name__, instance_relative_config=True )
app.logger.setLevel( logging.INFO )

@app.route( "/" )
def mainpage():
    return flask.render_template( 'fom.html' )


@app.route( "/calculate", methods=['POST'] )
def calculate( logger=logging.getLogger("main") ):
    response = {}
    result = subprocess.run( workdir / "snflux", capture_output=True )
    if result.returncode != 0:
        logger.error( f"snflux returned {result.returncode}\n"
                      f"STDOUT\n------\n{result.stdout}\n"
                      f"STDERR\n------\n{result.stderr}" )
        response = { 'error': 'Error return from snflux' };
        return response

    prism.prism( flask.request.json )
    sninvar = pathlib.Path( workdir / "sninvarim.txt" )
    sninvar.replace( workdir / "sninvar.txt" )
    
    result = subprocess.run( workdir / "snvar", capture_output=True )
    if result.returncode != 0:
        logger.error( f"snvar returned {result.returncode}\n"
                      f"STDOUT\n------\n{result.stdout}\n"
                      f"STDERR\n------\n{result.stderr}\n" )
        response = { 'error': 'Error return from snvar' }
        return response

    with open( workdir / "snoutput.txt", "r" ) as ifp:
        snoutput = ifp.read()

    return { 'status': 'ok', 'text': snoutput }


# def create_app():

# app.add_url_rule( "/", view_func=mainpage )
# app.add_url_rule( "/calculate", view_func=lambda: calculate( logger=app.logger ) )
    
# return app


