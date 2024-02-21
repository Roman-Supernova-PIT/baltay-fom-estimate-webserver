import { fom } from "./fom.js"

// **********************************************************************
// **********************************************************************
// **********************************************************************
// Here is the thing that will make the code run when the document has loaded

fom.started = false

// console.log("About to window.setInterval...");
fom.init_interval = window.setInterval(
    function()
    {
        var requestdata, renderer
        
        if (document.readyState == "complete")
        {
            // console.log( "document.readyState is complete" );
            if ( !fom.started )
            {
                fom.started = true;
                window.clearInterval( fom.init_interval );
                renderer = new fom.Context();
                renderer.renderpage();
            }
        }
    },
    100);

export { }
