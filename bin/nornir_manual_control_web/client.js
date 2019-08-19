/*jshint node:true, esnext: true, laxcomma: true, smarttabs: true */
'use strict';

var nanomsg = require('nanomsg')
  , http    = require('http')
  , fs      = require('fs')
  , path    = require('path')
  , socket  = nanomsg.socket('pair')
  , methods = {};

socket.bind('ws://*:3001');

function handler( req, res ){
    var method = req.method.toLowerCase();
    switch( method ){
        case 'get':
        case 'post':
            methods[ 'on_' + method ]( req, res ); 
            break;
        default:
            return _404( req, res );
    }
}

methods.on_post = function( req, res ){
    var data = "";

    req.on('data', function( buff ){
        data += buff.toString('utf8');
    });

    req.on('end', function(){
        socket.send( data );
        res.writeHead(201,{'Content-Type': 'application/json'});
        res.end( );
    });
};

methods.on_get = function( req, res ){
    res.writeHead(200,{'Content-Type':'text/html'});
    fs.createReadStream( path.join(__dirname, 'index.html')).pipe( res );
};

function _404( req, res ){
    res.writeHead(405, {'Content-Type':'text/plain'});
    res.end();
}

http.createServer( handler ).listen(3000);