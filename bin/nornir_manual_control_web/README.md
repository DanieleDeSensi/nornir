Based on https://github.com/denilsonsa/html5-knob, customized for the purpose

Steps to let it run:
    1. install nodejs
    2. install npm
    3. run 'npm install nanomsg'
    4. run 'node client.js'
    5. open index.html with a browser and play.

server.cpp file is to be used for debugging purposes (to check data sent by the web interface).

It is composed by three part:
    - A web interface (index.html)
    - A javascript service (client.js)
    - The nornir manager (to be run with strategySelection = STRATEGY_SELECTION_MANUAL_WEB)

The web interface shows the knobs. When the values of the knobs are modified, the values are sent as string to the javascript service by using
a POST request. The javascript service, when the post request is received, will send the string to the nornir manager through websockets
by using nanomsg.

The web interface and the javascript service can run on different machines. In this case, the web interface must point to
the ip address of the javascript service. To do that, you should replace "localhost" in index.html with the address of the machine
running the javascript service.
In principle even the javascript service and the nornir manager could run on two different machines. However this has never
been tested and we strongly encourage you to run both of them on the same machine.
