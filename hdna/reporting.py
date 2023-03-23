import math
import pandas as pd
import numpy as np 
import plotly.graph_objects as go
import plotly.express as px

def upperapprox(num):
    order = math.floor(math.log10(num))
    return math.ceil(num/(10**order))*10**order

def lowerapprox(num):
    order = math.floor(math.log10(num))
    return math.floor(num/(10**order))*10**order

def num_ticks(lower, upper, tick_mag):
    return math.ceil((upper - lower + 1) / tick_mag)

def valplot(data, name, writepath = None, theme = 'dark', log=False, corr=False):
    
    THM = themetemplates(theme, 'scatter')

    sdata = data.sort_values('experimental')
    # sdata.set_index('seq', inplace=True, drop=False)
    try: sdata.drop('Unnamed: 0', axis=1)
    except KeyError: pass 
    X = list(sdata['experimental'].astype(float))
    if corr == False:
        Y = list(sdata['computational'].astype(float))
    else: Y = list(sdata['corrected'].astype(float))
    S = list(sdata['sequences'])
    I = list(sdata['index'])

    ERS = list(sdata['var'])

    if log == True:
        X = [np.log10(x) for x in X]
        Y = [np.log10(y) for y in Y]
        lowbound = lowerapprox(min(Y))
        topbound = upperapprox(max(X))
        nticks = 20
        spacing = (topbound-lowbound)/(nticks/2 +.5)
        xtitle = ''#r'$log_\text{10}(k_{exp}(seq)$'
        ytitle = ''#r'$log_\text{10}(k_{mod}(seq)$'
    else: 
        lowbound = 0
        topbound = 8e6#upperapprox(max(X))
        spacing = 1e6
        xtitle = ''#r'$k_{exp}(seq)$'
        ytitle = ''#r'$k_{mod}(seq)$'

    XLINE = np.linspace(lowbound, topbound)

    trace1 = go.Scatter(
        x = X, # My list of values for 'x'
        y = Y, # My list of values for 'y'
        error_y=dict(
            type='data',
            array = ERS,
            visible=True,
            color = THM['colordots'],
            thickness=1.5,
            width=3),
        marker = dict(color = 'white'),
        mode = 'markers',
        name = '',
        customdata = [f'{i}_{s}' for i, s in zip(I, S)],
        hovertemplate="""emp = %{x:.3e}
                     <br>mod = %{y:.3e}
                     <br>seq:  %{customdata} </b>"""
    )
    trace2 = go.Scatter(
        x = XLINE,
        y = XLINE,
        mode = 'lines',
        marker = dict(color = THM['colorline']),
        line = dict(dash = 'dash'),
        name = 'bisector'
    )
    layout = go.Layout(
        template=THM['template'],
        # title = f'Scatterplot for {name}',
        xaxis_title=xtitle,
        yaxis_title=ytitle,
        showlegend=False,
        autosize = False,
        width = 600,
        height = 600,
        margin = dict(
            b = 50,
            t = 50,
            l = 50,
            r = 50,
            pad = 0
        ),
        xaxis = dict(
            tickmode = 'array',
            tickvals = np.linspace(lowbound, topbound, num_ticks(lowbound, topbound, spacing)),
            showgrid = True
            ),
        yaxis = dict(
            tickmode = 'array',
            tickvals = np.linspace(lowbound, topbound, num_ticks(lowbound, topbound, spacing)),
            showgrid = True
        )
    )
    dados = [trace1, trace2]
    fig = go.Figure(data = dados, layout = layout)
    
    fig.update_xaxes(exponentformat="e", titlefont={'size': 22})
    fig.update_yaxes(exponentformat="e", titlefont={'size': 22})
    
    fig.update_layout(xaxis_range=[lowbound,topbound])
    fig.update_layout(yaxis_range=[lowbound,topbound])

    fig.update_traces(marker = dict(size=14,
                                    line=dict(width=2,
                                              color=THM['colordots'])),
                        selector=dict(mode='markers'))

    if writepath == None:
        fig.show()
        return fig 
    else:
        PATH = f'{writepath}/{name}'
        fig.write_html(f"{PATH}.html")




def histotime(data, fit_gamma, fit_exp, runtime, exp=None, mod=None, seq=None, nbins=150, name='timehist', writepath=None, theme='light'):
    
    if theme == 'dark':
        THM = {'template': 'plotly_dark',
                    'colorfit': 'coral',
                    'colorbin': 'white',
                    'fitwidth': 3}
    if theme == 'light':
        THM = {'template': 'ggplot2',
                    'colorfit': 'royalblue',
                    'colorbin': 'tomato',
                    'fitwidth': 2.5}
    
    nbins = nbins
    X = np.linspace(0,runtime,500)
    data1 = go.Scatter(
                x=X, 
                y=fit_gamma.pdf(X),
                name='beta distribution fit',
                line = dict(
                    color = THM['colorfit'],
                    width =5))#THM['fitwidth']))
    data3 = go.Scatter(
                x=X, 
                y=fit_exp.pdf(X),
                name='exponential distribution fit',
                line = dict(
                    dash = 'dot',
                    color = '#463F3A',
                    width = THM['fitwidth']))
    data2 = go.Histogram(
                x=data, 
                histnorm='probability density', 
                name='Simulation FPT distribution',
                nbinsx=nbins,
                marker = dict(color = THM['colorbin']))
    
    layout = go.Layout(
                # title="""seq %-s<br>mod: %-.2e<br>exp: %-.2e""" % (name.split('_')[0], float(mod), float(exp)),
                width=800,
                height=800,
                xaxis_title="simulation time",
                yaxis_title="probability density",
                template='presentation',#THM['template'],
                xaxis=dict(tickformat=".1e"),
                yaxis=dict(tickformat=".1e"))

    fig = go.Figure(data=[data1, data2, data3], layout=layout)
    fig.update_traces(opacity=0.95)
    fig.update_layout(legend=dict(
        orientation='h',
        yanchor='bottom',
        y=1.03,
        xanchor='right',
        x = 1
    ))

    if writepath == None:
        fig.show()
    else:
        PATH = f'{writepath}/{name}'
        fig.write_html(f"{PATH}.html")
    


def themetemplates(choice, kind):
    if kind == 'scatter' or 'percent':
        if choice == 'dark':
            return {'template': 'plotly_dark',
                    'colordots': 'tan',
                    'colorline': 'coral',
                    'linewidth': 2}
        if choice == 'light':
            return {'template': 'presentation',
                    'colordots': 'tomato',
                    'colorline': 'royalblue',
                    'linewidth': 2}

        

def percomplot(fpts, writepath=None, name=None, theme='dark'):

    THM = themetemplates(theme, 'percent')

    fig = px.ecdf(fpts, ecdfnorm='percent')
    fig.update_layout(
        width = 800,
        height = 500,
        title="percent completion by time",
        xaxis_title="simulation time",
        yaxis_title="percent",
        yaxis_ticksuffix = "%",
        xaxis=dict(tickformat=".1e"),
        showlegend=False,
        template=THM['template']
        )
    
    fig.data[0].line.color = THM['colorline']
    fig.data[0].line.width = THM['linewidth']

    if writepath == None:
        fig.show()
    else:
        PATH = f'{writepath}/{name}'
        fig.write_html(f"{PATH}.html")

class Tee(object):
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush() # If you want the output to be visible immediately
    def flush(self):
        for f in self.files:
            f.flush()
    def close(self):
        for f in self.files:
            f.close()







    #         def violinplot(data, name, writepath = None, theme = 'dark', log=False, corr=False):
    
    # THM = themetemplates(theme, 'scatter')

    # sdata = data.sort_values('experimental')
    # # sdata.set_index('seq', inplace=True, drop=False)
    # try: sdata.drop('Unnamed: 0', axis=1)
    # except KeyError: pass 
    # X = list(sdata['experimental'].astype(float))
    # if corr == False:
    #     Y = list(sdata['computational'].astype(float))
    # else: Y = list(sdata['corrected'].astype(float))
    # S = list(sdata['sequences'])
    # I = list(sdata['index'])

    # ERS = list(data['var'])

    # if log == True:
    #     X = [np.log10(x) for x in X]
    #     Y = [np.log10(y) for y in Y]
    #     lowbound = lowerapprox(min(Y))
    #     topbound = upperapprox(max(Y))
    #     nticks = 20
    #     spacing = (topbound-lowbound)/(nticks/2 +.5)
    #     xtitle = ''#r'$log_\text{10}(k_{exp}(seq)$'
    #     ytitle = ''#r'$log_\text{10}(k_{mod}(seq)$'
    # else: 
    #     lowbound = 0
    #     topbound = upperapprox(max(Y))
    #     spacing = 1e6
    #     xtitle = ''#r'$k_{exp}(seq)$'
    #     ytitle = ''#r'$k_{mod}(seq)$'

    # XLINE = np.linspace(lowbound, topbound)

    # trace1 = go.Scatter(
    #     x = X, # My list of values for 'x'
    #     y = Y, # My list of values for 'y'
    #     error_y=dict(
    #         type='data',
    #         array = ERS,
    #         visible=True,
    #         color = 'green',
    #         thickness=1.5,
    #         width=3),
    #     marker = dict(color = THM['colordots']),
    #     mode = 'markers',
    #     name = '',
    #     customdata = [f'{i}_{s}' for i, s in zip(I, S)],
    #     hovertemplate="""emp = %{x:.3e}
    #                  <br>mod = %{y:.3e}
    #                  <br>seq:  %{customdata} </b>"""

    # )
    # trace2 = go.Scatter(
    #     x = XLINE,
    #     y = XLINE,
    #     mode = 'lines',
    #     marker = dict(color = THM['colorline']),
    #     line = dict(dash = 'dash'),
    #     name = 'bisector'
    # )
    # layout = go.Layout(
    #     template=THM['template'],
    #     title = f'Scatterplot for {name}',
    #     xaxis_title=xtitle,
    #     yaxis_title=ytitle,
    #     showlegend=False,
    #     autosize = False,
    #     width = 600,
    #     height = 600,
    #     margin = dict(
    #         b = 50,
    #         t = 50,
    #         l = 50,
    #         r = 50,
    #         pad = 0
    #     ),
    #     xaxis = dict(
    #         tickmode = 'array',
    #         tickvals = np.linspace(lowbound, topbound, num_ticks(lowbound, topbound, spacing)),
    #         showgrid = True
    #         ),
    #     yaxis = dict(
    #         tickmode = 'array',
    #         tickvals = np.linspace(lowbound, topbound, num_ticks(lowbound, topbound, spacing)),
    #         showgrid = True
    #     )
    # )
    # dados = [trace1, trace2]
    # fig = go.Figure(data = dados, layout = layout)
    
    # fig.update_xaxes(exponentformat="e", titlefont={'size': 22})
    # fig.update_yaxes(exponentformat="e", titlefont={'size': 22})
    
    # fig.update_layout(xaxis_range=[lowbound,topbound])
    # fig.update_layout(yaxis_range=[lowbound,topbound])

    # fig.update_traces(marker = dict(size=12,
    #                                 line=dict(width=2,
    #                                           color='DarkSlateGrey')),
    #                     selector=dict(mode='markers'))

    # if writepath == None:
    #     fig.show()
    #     return fig 
    # else:
    #     PATH = f'{writepath}/{name}'
    #     fig.write_html(f"{PATH}.html")