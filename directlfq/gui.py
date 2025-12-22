import sys
import os
import time
import logging
import pandas as pd

# visualization
import panel as pn
import bokeh.server.views.ws
import directlfq.dashboard_parts as dashboard_parts


def get_css_style(
    file_name="dashboard_style.css",
    directory=os.path.join(
        os.path.dirname(__file__),
        "style"
    )
):
    file = os.path.join(
        directory,
        file_name
    )
    with open(file) as f:
        return f.read()


def init_panel():
    pn.extension(raw_css=[get_css_style()])
    pn.extension('plotly')


# style
init_panel()


class GUI(object):
    # TODO: import from alphabase

    def __init__(
        self,
        name,
        github_url,
        run_in_background=False,
        automatic_close=True,
        port=None,
    ):
        self.name = name
        self.tab_counter = 0
        self.header = dashboard_parts.HeaderWidget(
            name,
            os.path.join(
                os.path.dirname(__file__),
                "img",
            ),
            github_url
        )
        self.layout = pn.Column(
            self.header.create(),
            sizing_mode='stretch_width',
            min_width=1270
        )
        self.run_in_background = run_in_background
        self.automatic_close = automatic_close
        self.port = int(port) if port else port

    def start_server(self, run_in_background=False):
        if self.automatic_close:
            bokeh_ws_handler = bokeh.server.views.ws.WSHandler
            self.bokeh_server_open = bokeh_ws_handler.open
            bokeh_ws_handler.open = self.__open_browser_tab(
                self.bokeh_server_open
            )
            self.bokeh_server_on_close = bokeh_ws_handler.on_close
            bokeh_ws_handler.on_close = self.__close_browser_tab(
                self.bokeh_server_on_close
            )

        port_param = {} if self.port is None else {"port": self.port}
        self.server = self.layout.show(threaded=True, title=self.name, **port_param)
        if not run_in_background:
            self.server.join()
        elif not self.run_in_background:
            self.server.join()

    def __open_browser_tab(self, func):
        def wrapper(*args, **kwargs):
            self.tab_counter += 1
            return func(*args, **kwargs)
        return wrapper

    def __close_browser_tab(self, func):
        def wrapper(*args, **kwargs):
            self.tab_counter -= 1
            return_value = func(*args, **kwargs)
            if self.tab_counter == 0:
                self.stop_server()
            return return_value
        return wrapper

    def stop_server(self):
        logging.info("Stopping server...")
        self.server.stop()
        if self.automatic_close:
            bokeh_ws_handler = bokeh.server.views.ws.WSHandler
            bokeh_ws_handler.open = self.bokeh_server_open
            bokeh_ws_handler.on_close = self.bokeh_server_on_close


class AlphaQuantGUI(GUI):
    # TODO: docstring
    def __init__(self, start_server=False, port=None):
        super().__init__(
            name="directLFQ",
            github_url='https://github.com/MannLabs/directLFQ',
            port=port,
        )
        self.project_description = """### directLFQ provides ratio-based normalization and protein intensity estimation for small and very large numbers of proteomes."""
        self.manual_path = os.path.join(
            os.path.dirname(__file__),
            "docs",
            'Empty_manual.pdf'
        )
        self.main_widget = dashboard_parts.MainWidget(
            self.project_description,
           # self.manual_path
        )

        # ERROR/WARNING MESSAGES
        self.error_message_upload = "The selected file can't be uploaded. Please check the instructions for data uploading."

        self.run_pipeline = dashboard_parts.RunPipeline()
        #self.tabs = dashboard_parts.Tabs(self.run_pipeline)

        self.layout += [
            self.main_widget.create(),
            self.run_pipeline.create(),
         #   self.tabs.create(),
        ]
        if start_server:
            self.start_server()


def run(port=None):
    AlphaQuantGUI(start_server=True, port=port)


if __name__ == '__main__':
    run()
