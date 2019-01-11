"""
This module contains the methods used to store the results.

3 different types of output are available:

    - **tsv** file
    - **png** graph: uses the *tsv* file and matplotlib
    - **html** graph: uses the *tsv* file and bokeh
"""

import gzip
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging
import pandas as pd
import numpy as np
from bokeh import __version__ as bokeh_version
from bokeh.plotting import output_notebook
from bokeh.embed import components
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import HoverTool, PrintfTickFormatter, CustomJS, Circle
from bokeh.models.widgets.inputs import TextInput
from bokeh.layouts import column, widgetbox

from simReg import __logger_name__

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

logger = logging.getLogger(__logger_name__)



class QQPlot(object):
    """
    QQPlot

    Args:
        input_file: tsv file with the data
        cutoff (bool): add cutoffs to the figure
        rename_fields (dict): column names from the input file can be renamed providing a dictionary {old_name : new_name}
        extra_fields (list): list of column names that want to be passed to the figure data. Need for example to
            search by them.
    """

    def __init__(self, input_file, cutoff=True, rename_fields=None, extra_fields=None):
        basic_tools = "pan,box_zoom,resize,wheel_zoom,reset,previewsave,crosshair"
        self.figure = figure(width=600, plot_height=600, tools=basic_tools, toolbar_location="above")
        # labels
        self.figure.xaxis.axis_label = 'Expected p-values'
        self.figure.yaxis.axis_label = 'Observed p-values'
        self.figure.axis.major_label_text_font_size = '14pt'
        #self.figure.xaxis[0].formatter = PrintfTickFormatter(format="%.1f") #%.0e
        #self.figure.yaxis[0].formatter = PrintfTickFormatter(format="%.1f")  # 1e-%1.4f
        #self.figure.xaxis[0].major_label_orientation = pi / 4

        self.layout = None

        self.__load_values(input_file, rename_fields=rename_fields)

        self.__create_basic_plot(extra_fields=extra_fields, cutoff=cutoff)

    def __load_values(self, input_file, rename_fields=None):
        """
        Load the values for the plot as a :obj:`~pandas.DataFrame`

        Args:
            input_file: tsv file with data
            rename_fields (dict): pair of old column name - new column name

        """

        self.data = pd.read_csv(input_file, header=0, sep="\t")
        if rename_fields is not None:
            self.data = self.data.rename(columns = rename_fields)

    def __create_basic_plot(self, extra_fields=None, cutoff = True):
        """
        Creates the basic plot.

        Args:
            extra_fields (list): list of extra fields containen by each point
            cutoff (bool): add cutoff

        """

        # Default settings
        min_samples = 2
        min_pvalue = 10000

        colors = ['royalblue', 'blue']

        self.data['observed'] = self.data['pvalue'].map(lambda x: -np.log10(x) if x > 0 else -np.log10(1 / min_pvalue))
        self.data['color'] = self.data['num_samples'].map(lambda x: colors[1] if x >= min_samples else colors[0])
        self.data['alpha'] = self.data['num_samples'].map(lambda x: 0.7 if x >= min_samples else 0.3)

        self.data.sort_values(by=['observed'], inplace=True)
        exp_pvalues = -np.log10(np.arange(1, len(self.data) + 1) / float(len(self.data)))
        exp_pvalues.sort()
        self.data['expected'] = exp_pvalues

        dict_for_source = dict(
            x=self.data['expected'].tolist(),
            y=self.data['observed'].tolist(),
            color=self.data['color'].tolist(),
            alpha=self.data['alpha'].tolist(),
            pvalue=[str(x) for x in self.data["pvalue"]],
            qvalue=[str(x) for x in self.data["qvalue"]]
        )
        if extra_fields is not None:
            for field in extra_fields:
                dict_for_source[field] = self.data[field].tolist()

        self._source = ColumnDataSource(data=dict_for_source)

        # Plot the first set of data
        if len(self.data['expected']) > 0:
            invisible_circle = Circle(x='x', y='y', fill_color='color', fill_alpha='alpha', line_color=None, size=10)
            visible_circle = Circle(x='x', y='y', fill_color='color', fill_alpha=0.9, line_color='red', size=10)
            self.glyph = self.figure.add_glyph(self._source, invisible_circle, selection_glyph=visible_circle,
                                           nonselection_glyph=invisible_circle)

        # Get the maximum pvalues (observed and expected)
        max_x = float(self.data[['expected']].apply(np.max))
        max_y = float(self.data[['observed']].apply(np.max))

        # Give some extra space (+-10%)
        max_x *= 1.1
        # Give some extra space (+-20%)
        max_y *= 1.2

        # Add a dashed diagonal from (min_x, max_x) to (min_y, max_y)
        self.figure.line(np.linspace(0, min(max_x, max_y)),
                       np.linspace(0, min(max_x, max_y)),
                       color='red', line_width=2, line_dash=[5, 5])

        # Set the grid
        self.figure.grid.grid_line_alpha = 0.8
        self.figure.grid.grid_line_dash = [6, 4]

        #cutoff
        if cutoff:
            # FDR
            for fdr_cutoff, fdr_color in zip((0.25, 0.1), ('green', 'red')):
                fdr = self.data[self.data['qvalue'] < fdr_cutoff]['observed']
                if len(fdr) > 0:
                    fdr_y = np.min(fdr)
                    fdr_x = np.min(self.data[self.data['observed'] == fdr_y]['expected'])
                    self.figure.line((fdr_x - max_x * 0.025, fdr_x + max_x * 0.025), (fdr_y, fdr_y),
                                   color=fdr_color, line_width=2)

    def add_tooltip(self):
        """
        Adds tooltip to show the parameters of each glyph in the figure
        """
        tooltip = \
            """
            <div>
                <span style='font-size: 17px; font-weight: bold;'>@HugoID</span>
                <span style='font-size: 15px; color: #966;'>[@EnsemblID]</span>
            </div>
            <div>
                <span style='font-size: 15px;'>p, q-value</span>
                <span style='font-size: 10px; color: #696;'>(@pvalue, @qvalue)</span>
            </div>
            </br>

            """
        self.figure.add_tools(HoverTool(tooltips=tooltip, mode='mouse', renderers=[self.glyph]))

    def add_tooltip_enhanced(self):
        """
        The tooltip is shown via JavaScript to avoid been block in areas with a
        high density of points
        """
        tooltip = \
            """
            "<div> \\
                <span style='font-size: 17px; font-weight: bold;'>\"+s.HugoID[index]+\"</span> \\
                <span style='font-size: 15px; color: #966;'>[\"+s.EnsemblID[index]+\"]</span> \\
            </div> \\
            <div> \\
                <span style='font-size: 15px;'>p, q-value</span> \\
                <span style='font-size: 10px; color: #696;'>(\"+s.pvalue[index]+\", \"+s.qvalue[index]+\")</span> \\
            </div> \\
            </br>"
            """

        code = \
            """
            if (cb_data.index['1d'].indices.length > 0) {
                var s = source.get('data');
                if ( $( ".bk-tooltip.bk-tooltip-custom.bk-left" ).length == 0 ){
                    $( ".bk-canvas-overlays" ).append( '<div class="bk-tooltip bk-tooltip-custom bk-left" style="z-index: 1010; top: 0px; left: 0px; display: block;"></div>' );
                }
                var inner = "<div><div><div>";
                for(i in cb_data.index['1d'].indices){
                    index = cb_data.index['1d'].indices[i];
                    if (i > 2) break;
                     inner = inner + """ + tooltip + """
                }
                if (i>2) {
                    inner = inner + "<div><span style='font-size: 15px; color: #009;'>TOTAL OF: " + cb_data.index['1d'].indices.length + "</span></div>";
                }
                inner = inner + "</div></div></div>";
                $('.bk-tooltip.bk-tooltip-custom.bk-left')[0].innerHTML = inner;
                $('.bk-tooltip.bk-tooltip-custom.bk-left').attr('style', 'left:' + (cb_data.geometry.sx+10) + 'px; top:' + (cb_data.geometry.sy-5-$('.bk-tooltip.bk-tooltip-custom.bk-left').height()/2) + 'px; z-index: 1010; display: block;');
            }else {
                $("div").remove(".bk-tooltip.bk-tooltip-custom.bk-left");
            }

            """

        callback = CustomJS(args={'source': self._source}, code=code)
        self.figure.add_tools(HoverTool(tooltips=None, callback=callback, renderers=[self.glyph], mode='mouse'))

    def add_search_widget(self, fields):
        """
        Add text input for each field.

        Args:
            fields (:obj:`str` or :obj:`list`): list of fields to do a search.

        """
        if not fields:
            return

        if isinstance(fields, str):
            callback = CustomJS(args=dict(source=self._source),
                                code=self.__get_search_box_code(fields))
            widgets = TextInput(value="", name=fields, title=fields,
                                callback = callback)
        else:
            widgets = []
            for field in fields:
                callback = CustomJS(args=dict(source=self._source),
                                    code=self.__get_search_box_code(field))
                text_input = TextInput(value="", name=field, title=field,
                                       callback=callback)
                widgets.append(text_input)

        widgets = widgetbox(widgets)

        self.layout = column(widgets, self.figure)


    @staticmethod
    def __get_search_box_code(field):
        """
        Code for a search box widget

        Args:
            field (str): field name

        Returns:
            str: Javascript code

        """
        code = \
            """
            origSearch = cb_obj.get('value');
            search = origSearch.toUpperCase();
            var selected = source.get('selected')['1d'].indices;
            searcher = source.get('data').""" + field + """
            selected.length = 0
            for(index in searcher) {
                if ( searcher[index].toUpperCase().indexOf(search) > -1) {
                    selected.push(index);
                }
            }
            if (selected.length == 0) {
                swal("Error", "Value not found: '"+origSearch+"'", "error")
            }
            source.trigger('change');
            """

        return code

    def show(self, output_path, showit=True, notebook=False):
        """
        Show the figure

        Args:
            output_path: file where to store the figure
            showit (bool): the figure is displayed (widgets and the like are not shown) or is fully saved. Defaults to True.
            notebook (bool): if is is called form a notebook or not. Defaults to False.

        """
        # Import modules
        if notebook:
            output_notebook()
        else:
            if output_path is not None:
                output_file(output_path)

        if self.layout is not None:
            layout = self.layout
        else:
            layout = self.figure

        if showit:
            show(self.figure)
        else:
            script, div = components(layout)
            html = """  <!DOCTYPE html>
                        <html>
                        <script src="https://code.jquery.com/jquery-1.11.3.min.js">
                        </script><script src="https://cdnjs.cloudflare.com/ajax/libs/sweetalert/1.1.3/sweetalert.min.js"></script>
                        <link rel="stylesheet" type="text/css" href="https://cdnjs.cloudflare.com/ajax/libs/sweetalert/1.1.3/sweetalert.min.css">
                        <link href="http://cdn.pydata.org/bokeh/release/bokeh-""" + bokeh_version + """.min.css" rel="stylesheet" type="text/css">
                        <script src="http://cdn.pydata.org/bokeh/release/bokeh-""" + bokeh_version + """.min.js"></script>
                        <link href="http://cdn.pydata.org/bokeh/release/bokeh-widgets-""" + bokeh_version + """.min.css" rel="stylesheet" type="text/css">
                        <script src="http://cdn.pydata.org/bokeh/release/bokeh-widgets-""" + bokeh_version + """.min.js"></script>
                        <style>.bk-root .bk-toolbar-above {right: auto;}</style>\n""" + \
                   script + \
                   """<body>\n""" + \
                   div + """
                        </body>
                        </html>"""
            text_file = open(output_path, "w")
            text_file.write(html)
            text_file.close()


def eliminate_duplicates(df):

    colors = df['color'].tolist()
    x = df['x'].tolist()[0]
    y = df['y'].tolist()[0]
    if 'r' in colors:
        return x, y, 'red'
    elif 'g' in colors:
        return x, y, 'green'
    else:
        return x, y, 'grey'


def add_symbol(df):
    ensemble_file = os.path.join(__location__, "ensembl_genes_75.txt.gz")
    gene_conversion = {line.split("\t")[0]: line.strip().split("\t")[-1]
                       for line in gzip.open(ensemble_file, 'rt').readlines()}
    gene_symbols = df.GENE_ID.apply(lambda e: gene_conversion.get(e, e))
    df.SYMBOL.fillna(gene_symbols, inplace=True)
    return df


def store_png(input_file, output_file, showit=False):
    """
    Creates a figure from the resutls.

    Args:
        input_file: tsv file with the results
        output_file: file where to store the figure
        showit (bool): calls :func:`~matplotlib.pyplot.show` before returning.
            Defaults to False.

    """
    pvalue='P_VALUE'
    qvalue='Q_VALUE'
    min_samples = 2
    draw_greys = True
    annotate = True
    cut = True
    #############################

    MIN_PVALUE = 10000
    MAX_ANNOTATED_GENES = 50

    df = pd.read_csv(input_file, header=0, sep="\t")
    df.dropna(subset=[pvalue], inplace=True)

    # Define the shape of the figure
    NROW = 1
    NCOL = 1

    fig = plt.figure(figsize=(6, 6))
    axs = [plt.subplot2grid((NROW, NCOL), (item // NCOL, item % NCOL)) for item in range(NROW * NCOL)]

    # Plot is on the right or on the left?
    dx_side = True
    ax = axs[0]

    colors = ['royalblue', 'blue']
    obs_pvalues = df[pvalue].map(lambda x: -np.log10(x) if x > 0 else -np.log10(1 / MIN_PVALUE))
    #gens with num_samples >= min_samples have different values for color and alpha
    obs_color = colors[1]
    obs_alpha = 0.7

    data = pd.DataFrame({'HugoID': df['SYMBOL'],
                         'observed': obs_pvalues,
                         'color': obs_color,
                         'alpha': obs_alpha,
                         'fdr': df[pvalue] if pvalue is not None else 1
                         })

    data.sort_values(by=['observed'], inplace=True)
    exp_pvalues = -np.log10(np.arange(1, len(data) + 1) / float(len(data)))
    exp_pvalues.sort()
    data['expected'] = exp_pvalues

    # Get the maximum pvalues (observed and expected)
    max_x = float(data[['expected']].apply(np.max))
    max_y = float(data[['observed']].apply(np.max))

    # Give some extra space (+-5%)
    min_x = max_x * -0.05
    min_y = max_y * -0.05
    max_x *= 1.1
    max_y *= 1.1

    grey = data[data['color'] == colors[0]]
    blue = data[data['color'] == colors[1]]

    # Plot the data
    if draw_greys and len(grey['expected']) > 0:
        ax.scatter(grey['expected'].tolist(),
                   grey['observed'].tolist(),
                   color=grey['color'].tolist(),
                   alpha=grey['alpha'].tolist()[0],
                   s=30)

    if len(blue['expected']) > 0:
        ax.scatter(blue['expected'].tolist(),
                   blue['observed'].tolist(),
                   color=blue['color'].tolist(),
                   alpha=blue['alpha'].tolist()[0],
                   s=30)

    # Get the data that are significant with a FDR < 0.1 and FDR 0.25
    genes_to_annotate = []
    for fdr_cutoff, fdr_color in zip((0.05, 0.01), ('g-', 'r-')):
        fdr = data[data['fdr'] < fdr_cutoff]['observed']
        if len(fdr) > 0:
            fdr_y = np.min(fdr)
            fdr_x = np.min(data[data['observed'] == fdr_y]['expected'])
            ax.plot((fdr_x - max_x * 0.025, fdr_x + max_x * 0.025), (fdr_y, fdr_y), fdr_color)
            # Add the name of the significant genes
            genes = data[(data['observed'] >= fdr_y) & (data['expected'] >= fdr_x)]
            for count, line in genes.iterrows():
                if line['color'] == colors[0]:
                    continue
                genes_to_annotate.append({'x': line['expected'],
                                          'y': line['observed'],
                                          'HugoID': line['HugoID'],
                                          'color': colors[0] if line['color'] == colors[0] else fdr_color[0]})

    # Annotate the genes
    genes_annotated = 0
    if annotate and len(genes_to_annotate) > 0:
        genes_to_annotate = pd.DataFrame(genes_to_annotate)

        # Get rid of redundant genes
        grouped = genes_to_annotate.groupby(['HugoID'])
        grouped = grouped.apply(eliminate_duplicates)
        grouped = pd.DataFrame({'HugoID': grouped.index.tolist(),
                                'x': [x[0] for x in grouped],
                                'y': [y[1] for y in grouped],
                                'color': [c[2] for c in grouped]})
        grouped.sort_values(by=['y', 'x'], inplace=True, ascending=[False, False])

        x_text = max_x * 1.1 if dx_side is True else min_x - (max_x * 0.2)
        y_text = np.floor(max_y)
        distance_between_genes = max_y * 0.05
        for count, line in grouped.iterrows():
            x, y = line['x'], line['y']
            ax.annotate(line['HugoID'], xy=(x, y), xytext=(x_text, y_text),
                        arrowprops=dict(facecolor="black", shrink=0.03,
                                        width=1, headwidth=6, alpha=0.3),
                        horizontalalignment="left", verticalalignment="center",
                        color=line['color'],
                        weight = 'normal'
                        )
            y_text -= distance_between_genes

            # This avoid the crash for ValueError: width and height must each be below 32768
            genes_annotated += 1
            if genes_annotated >= MAX_ANNOTATED_GENES:
                logger.warning("Annotations cut to %d genes", MAX_ANNOTATED_GENES)
                break

    # Add labels
    ax.set_xlabel("expected pvalues")
    ax.set_ylabel("observed pvalues")

    # Add the dashed diagonal
    ax.plot(np.linspace(0, np.floor(max(max_x, max_y))),
            np.linspace(0, np.floor(max(max_x, max_y))),
            'r--')
    ax.grid(True)

    # Redefine the limits of the plot
    if cut:
        ax.set_xlim(min_x, max_x)
        ax.set_ylim(min_y, max_y)

    # Set the title: project, cancer_type
    #ax.set_title(title)

    # Adjust the plot
    try:
        plt.tight_layout()
    except ValueError as e:
        logger.warning('Ignoring tight_layout()')

    # Save the plot
    if output_file:
        plt.savefig(output_file, bbox_inches='tight')

    # Show the plot
    if showit:
        plt.show()

    # Close the figure
    plt.close()


def store_html(input_file, output_path):
    """
    Create the QQPlot and save it.

    Args:
        input_file: tsv filw with the data
        output_path: file where to store the graph
        showit (bool): defaults to False. See :meth:`~simReg.store.QQPlot.show`.

    """

    search_by_fields = ['HugoID', 'EnsemblID']

    qqp = QQPlot(input_file=input_file, rename_fields={'SAMPLES': 'num_samples','SYMBOL': 'HugoID', 'GENE_ID': 'EnsemblID', 'P_VALUE': 'pvalue', 'Q_VALUE': 'qvalue'}, extra_fields=search_by_fields, cutoff=True)

    qqp.add_tooltip_enhanced()

    qqp.add_search_widget(search_by_fields)

    qqp.show(output_path = output_path, notebook = False, showit=False)


def store_tsv(results, result_file):
    """
    Saves the results in a tsv file sorted by pvalue

    Args:
        results (:obj:`~pandas.DataFrame`): results of the analysis
        result_file: file where to store the results

    """
    results.index.names = ['GENE_ID']
    results.sort_values(by='pvalue', inplace=True)
    fields = ['muts', 'muts_recurrence', 'samples_mut', 'pvalue', 'qvalue', 'pvalue_neg', 'qvalue_neg', 'snps', 'mnps', 'indels', 'symbol']
    df = results[fields].copy()
    df.reset_index(inplace=True)
    df.rename(columns={'muts': 'MUTS', 'muts_recurrence': 'MUTS_RECURRENCE', 'samples_mut': 'SAMPLES',
                       'pvalue': 'P_VALUE', 'qvalue': 'Q_VALUE','pvalue_neg': 'P_VALUE_NEG',
                       'qvalue_neg': 'Q_VALUE_NEG', 'snps': 'SNP', 'mnps':'MNP', 'indels': 'INDELS',
                       'symbol': 'SYMBOL'}, inplace=True)
    df = add_symbol(df)

    with open(result_file, 'wt') as fd:
        df.to_csv(fd, sep="\t", header=True, index=False)
