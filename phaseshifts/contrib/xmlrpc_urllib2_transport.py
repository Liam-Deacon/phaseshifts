#!/bin/env python
"""urllib2-based transport class for xmlrpclib.py (with test code).

Written from scratch but inspired by xmlrpc_urllib_transport.py file from http://starship.python.net/crew/jjkunce/ by jjk.

A. Ellerton 2006-07-06

Testing with Python 2.4 on Windows and Linux, with/without a corporate proxy in place.

****************************
*** USE AT YOUR OWN RISK ***
****************************
"""

import xmlrpclib


class ProxyTransport(xmlrpclib.Transport):
    """Provides an XMl-RPC transport routing via a http proxy.
    This is done by using urllib2, which in turn uses the environment
    varable http_proxy and whatever else it is built to use (e.g. the
    windows registry).
    NOTE: the environment variable http_proxy should be set correctly.
    See checkProxySetting() below.
    Written from scratch but inspired by xmlrpc_urllib_transport.py
    file from http://starship.python.net/crew/jjkunce/ by jjk.
    A. Ellerton 2006-07-06"""

    # redefine parse_response and _parse_response(taken from python 2.6)
    # to work with python 2.7

    def parse_response(self, file):
        # compatibility interface
        return self._parse_response(file, None)

    def _parse_response(self, file, sock):
        # read response from input file/socket, and parse it

        p, u = self.getparser()

        while 1:
            if sock:
                response = sock.recv(1024)
            else:
                response = file.read(1024)
            if not response:
                break
            if self.verbose:
                print("body:", repr(response))
            p.feed(response)

        file.close()
        p.close()

        return u.close()

    def request(self, host, handler, request_body, verbose):
        import urllib2

        self.verbose = verbose
        url = "http://" + host + handler
        if self.verbose:
            "ProxyTransport URL: [%s]" % url

        request = urllib2.Request(url)
        request.add_data(request_body)
        # Note: 'Host' and 'Content-Length' are added automatically
        request.add_header("User-Agent", self.user_agent)
        request.add_header("Content-Type", "text/xml")  # Important

        proxy_handler = urllib2.ProxyHandler()
        opener = urllib2.build_opener(proxy_handler)
        f = opener.open(request)
        return self.parse_response(f)


def checkProxySetting():
    """If the variable 'http_proxy' is set, it will most likely be in one
    of these forms (not real host/ports):
    proxyhost:8080
    http://proxyhost:8080
    urlllib2 seems to require it to have 'http;//" at the start.
    This routine does that, and returns the transport for xmlrpc."""
    import os, re

    try:
        http_proxy = os.environ["http_proxy"]
    except KeyError:
        return

    # ensure the proxy has the 'http://' at the start
    #

    match = re.search(
        "(http://)?([\w/-/.]+):([\w/-/.]+)(\@)?([\w/-/.]+)?:?([\w/-/.]+)?", http_proxy
    )
    if not match:
        raise Exception("Proxy format not recognised: [%s]" % http_proxy)
    else:
        groups = match.groups()
        if not groups[3]:
            # proxy without authorization
            os.environ["http_proxy"] = "http://%s:%s" % (groups[1], groups(2))
        else:
            os.environ["http_proxy"] = "http://%s:%s@%s:%s" % (
                groups[1],
                groups[2],
                groups[4],
                groups[5],
            )

    return


def test():
    import sys, os

    def nextArg():
        try:
            return sys.argv.pop(1)
        except:
            return None

    checkProxySetting()

    url = nextArg() or "http://betty.userland.com"
    api = nextArg() or "examples.getStateName(32)"  # "examples.getStateList([1,2])"
    try:
        server = xmlrpclib.Server(url, transport=ProxyTransport())
        print("Url: %s" % url)

        try:
            print("Proxy: %s" % os.environ["http_proxy"])
        except KeyError:
            print("Proxy: (Apparently none)")

        print("API: %s" % api)
        r = eval("server.%s" % api)
        print("Result: ", r)

    except xmlrpclib.ProtocolError as err:
        print("Connection error: %s" % err)
    except xmlrpclib.Fault as err:
        print("Error: %s" % err)


if __name__ == "__main__":
    # run with no parameters for basic test case.
    test()
