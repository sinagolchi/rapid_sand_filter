import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import LineString

st.title('Biofilter Calculations')
st.caption('Developed by Sina Golchi for UW NSERC Chair in Water Treatment, 2022')

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(['Setup', 'Clean bed head loss' ,'Fluidization Q', 'Bed expansion by Q', 'Q for bed expansion', 'Pump sizing'])

with tab1:
    st.subheader('Sand bed configuration')
    c1,c2,c3 = st.columns(3)
    area = c1.number_input(label='Base area in [m^2]', value= 0.0079, step=0.0001,format="%.5f")
    total_length = c2.number_input(label='Total length of bed [m] (overflow position)', value=2.5)
    media_l = c2.number_input(label='Media depth at rest [m]', value=1.2)
    media_e = c3.number_input(label='Media porosity (0-1)', value= 0.40)
    media_rho = c3.number_input(label='Media density [kg/m^3]', value=2650)
    Kv = c2.number_input(label='Kv value for Ergun equation', value=115)
    Ki = c2.number_input(label='KI value for Ergun equation', value=2.5)
    d = c3.number_input(label='Effective size of media [mm]', value=0.6)
    temp = c1.number_input(label='Temperature [degree C]', value = 15)
    visc =  c1.number_input(label='Water dynamic viscosity [kg/m.s]',value=1.14E-3,format="%.5f")
    w_rho = c1.number_input(label='Water density [kg/m^3]', value=999)


    g = 9.81

    media_volume = area*media_l
    weight_sand = (1-media_e)*media_volume*media_rho
    st.metric(label='Media volume in M^3', value = np.around(media_volume,3))
    st.metric(label='Mass of sand in KG', value = np.around(weight_sand,3))

    st.markdown("""---""")
    st.subheader('Gravel bed configuration')
    grav_bed = st.checkbox(label='Account for head loss in gravel bed?')
    if grav_bed:
        c41, c42 = st.columns(2)
        grav_media_l = c41.number_input(label='Gravel depth[m]', value=0.3)
        grav_media_e = c42.number_input(label='Gravel porosity (0-1)', value=0.45)
        grav_Kv = c41.number_input(label='Gravel Kv value for Ergun equation', value=115)
        grav_Ki = c42.number_input(label='Gravel KI value for Ergun equation', value=2.5)
        grav_d = c41.number_input(label='Effective size of Gravel [mm]', value=3)

with tab2:
    c11,c12 = st.columns(2)
    Q_set_lpm = c11.number_input(label='Flow rate through bed [L/min]', value=10)
    v = Q_set_lpm*1.66667e-5/area
    c12.metric('Filtration rate [m/s]',value=np.around(v,5))
    def clean_bed_hlcalc(Kv,Ki,media_e,media_l,d,v):
        viscous_loss = ((Kv*((1-media_e)**2))*visc*media_l*v)/((media_e**3) * w_rho * g * (d * 1e-3) ** 2)
        inertia_loss = (Ki*(1-media_e)*(media_l)*v**2)/(media_e**3*g*(d*1e-3))
        return [viscous_loss, inertia_loss]

    visc_loss, inert_loss = clean_bed_hlcalc(Kv,Ki,media_e,media_l,d,v)
    c11.metric('Viscous head loss [m]', value=np.around(visc_loss,4))
    c12.metric('Inertia head loss [m]', value=np.around(inert_loss,4))
    st.metric('Total clean bed head loss [m]', value=np.around(visc_loss+inert_loss, 4))

    st.markdown("""___""")

    if grav_bed:
        grav_v_loss , grav_i_loss = clean_bed_hlcalc(grav_Kv,grav_Ki,grav_media_e,grav_media_l,grav_d,v)
        c11, c12,c13 = st.columns(3)
        c11.metric('Viscous head loss through gravel', value=str(np.around(grav_v_loss, 4)) + ' m')
        c12.metric('Inertia head loss through gravel', value=str(np.around(grav_i_loss, 4)) + ' m')
        c13.metric('Total clean bed head loss of gravel', value=str(np.around(grav_v_loss + grav_i_loss, 4)) + ' m')
        st.metric('Total bed head loss (Sand + Gravel)',value=str(np.around(visc_loss+inert_loss + grav_v_loss + grav_i_loss, 4)) + ' m')

    with st.expander('How does it work?'):
        st.markdown(r"""The clean bed head loss through a granular filter can be calculated by the Ergun equation (Ergun 1952):""")
        st.latex( r'h_L = K_V  \frac{(1-\epsilon)^2}{\epsilon^3} \frac{\mu L v}{\rho_W g d^2} + K_I \frac{1-\epsilon}{\epsilon^3}\frac{L v^2}{gd}')

with tab3:
    st.header('Fluidization head loss')
    fl_weight = (media_rho - w_rho)*(1-media_e)*area*media_l
    fl_force = fl_weight*g
    fl_hl = fl_force/(area*w_rho*g)
    c21, c22, c23 = st.columns(3)
    c21.metric('Fluidized weight of media [Kg]', value = np.around(fl_weight,4))
    c22.metric('Fluidized downward force of media [N]', value = np.around(fl_force,4))
    c23.metric('Fluidized head loss of media [m]', value = np.around(fl_hl,4))

    st.header('Minimum flow rate for fluidization')
    c31, c32 = st.columns(2)

    def back_wash_flow(Kv,Ki,d,media_e,media_rho):
        Beta = ((g*w_rho*(media_rho-w_rho)) * ((d * 1e-3)**3) * (media_e**3))/visc**2

        fl_re = ((-Kv*(1-media_e)) + np.sqrt((Kv**2 * (1-media_e)**2) + 4*Ki*Beta))/(2*Ki)

        fl_v = (visc*fl_re)/(w_rho*(d * 1e-3))

        fl_flow = (fl_v*area)/1.66667e-5
        return [fl_flow, fl_v, fl_re, Beta]


    fl_flow, fl_v, fl_re, Beta = back_wash_flow(Kv, Ki, d, media_e, media_rho)

    c31.metric('Backwash fluidization calculation factor', np.around(Beta, 4))
    c32.metric('Minimum Reynolds number for fluidization', np.around(fl_re, 4))
    c31.metric('Minimum filtration rate for fluidization (m/s)', np.around(fl_v, 4))
    c32.metric('Minimum flow rate for fluidization (L/min)', np.around(fl_flow, 4))

with tab4:
    st.subheader('How much the bed expand by a certain flow rate')
    c11,c12 = st.columns(2)
    Q_fl_lpm = c11.number_input(label='Flow rate through bed for fluidization [L/min]', value=10.0, step=0.1)
    v_fl = Q_fl_lpm*1.66667e-5/area
    c12.metric('Filtration rate [m/s]',value=np.around(v_fl,5))

    def bed_expansion_by_Q(v,kv,ki):
        X = ((visc*v)/(2*g*(media_rho-w_rho)*(d * 1e-3) ** 2))*(kv+((ki*w_rho*v*d * 1e-3)/visc))
        Y = (kv*visc*v)/(3*g*(media_rho-w_rho)*(d * 1e-3) ** 2)
        expand_e = np.cbrt(X+np.sqrt(X**2 + Y**3)) + np.cbrt(X-np.sqrt(X**2 + Y**3))
        expanded_depth = media_l*((1-media_e)/(1-expand_e))
        return [expanded_depth, expand_e]

    expanded_depth , expanded_e = bed_expansion_by_Q(v_fl,Kv,Ki)

    c11,c12 = st.columns(2)

    c11.metric('Media depth at rest', value=str(np.around(media_l, 3)) + ' m')
    c12.metric('Expanded depth',value=str(np.around(expanded_depth,3)) + ' m')
    c11.metric('Bed porosity at rest',value=str(np.around(media_e, 3)))
    c12.metric('Expanded bed porosity', value=str(np.around(expanded_e, 3)))
    st.metric('Total expansion', value=str(np.around((expanded_depth-media_l)/media_l *100, 3)) + ' %')
    if expanded_depth > total_length:
        st.warning(F'The backwash flow of {Q_fl_lpm} L/min will cause the bed to overflow! ')

with tab5:
    st.subheader('Calculating backwash flow for certain bed expansion')
    re_expand = st.number_input('Desired bed expansion [%]', min_value=10 ,value= 30, max_value=200)
    c41, c42 = st.columns(2)
    re_depth = media_l + (re_expand/100)*media_l
    re_e = 1 - ((media_l/re_depth)*(1-media_e))
    c41.metric('Media depth at rest', value=str(np.around(media_l, 3)) + ' m')
    c42.metric('Expanded media depth', value=str(np.around(re_depth, 3)) + ' m')
    c41.metric('Bed porosity at rest',value=str(np.around(media_e, 3)) )
    c42.metric('Expanded bed porosity', value=str(np.around(re_e, 3)) )

    res = back_wash_flow(Kv, Ki, d, re_e, media_rho)

    c31, c32 = st.columns(2)
    c31.metric(F'Filtration rate to achieve {str(re_expand)}% expansion',str(np.around(res[1], 4)) + ' m/s')
    c32.metric(F'Flow rate to achieve {str(re_expand)}% expansion', str(np.around(res[0], 4)) + ' L/min')

with tab6:
    st.subheader('Static and minor and major head loss')

    st.markdown('#### Minor and major losses')
    elev_hl = st.number_input('Elevation head loss (elevation to bed, not to top of column)', value=1.0, step=0.1)
    with st.expander('Pipe friction head loss configuration'):
        c51, c52, c53 = st.columns(3)
        pipe_diam = c51.number_input('Pipe diameter [mm]', value=20.0, step=0.1)
        pipe_eps = c52.number_input('Pipe roughness [mm]', value=0.003334, step=0.000001, format='%.6f')
        pipe_length = c53.number_input('Pipe length [m]', value=2.0, step=0.1)
        use_darcy = st.checkbox('Use Darcy-Weisbach friction factor')


    with st.expander('Fittings'):
        fitting_kvalue = st.number_input('Total K values of fittings',value=6.75)

    st.markdown('### Valve control (Extra head loss introduced to control the flow)')

    use_valve = st.checkbox('Use a valve to control the flow?')
    if use_valve:
        max_kv = st.number_input('Valve Kv value when fully open',value=1.47 , step=0.01)
        valve_Kv = st.select_slider('Valve Kv Value', options= np.linspace(max_kv,0.1,20))
    else:
        valve_Kv = None

    def system_pressure(flow,total_k,valve_Kv,use_valve,use_darcy):
        v_bed = flow * 1.66667E-5 / area
        v_pipe = flow * 1.66667E-5/ (3.14*((pipe_diam/1000)/2)**2)

        if use_darcy:
            reynolds_pipe = (w_rho*v_pipe*(pipe_diam/1000))/visc
            print(reynolds_pipe)
            s = 0.12363*reynolds_pipe*(pipe_eps/pipe_diam) + np.log(0.3984*reynolds_pipe)
            f = np.sqrt(1/( 0.8686*np.log((0.3984*reynolds_pipe)/(0.8686*s)**((s-0.645)/s+0.39))))
            k_pipe = f*(pipe_length/(pipe_diam/1000))
            total_k += k_pipe


        h_total = 0

        if grav_bed:
            h_total += (grav_v_loss +  grav_i_loss)

        if use_valve:
            flow_for_Kv = flow * 0.06 #m3/h
            delta_P_bar = (flow_for_Kv**2)/(valve_Kv**2)
            hl_valve = delta_P_bar*10.21
            h_total += hl_valve

        hl1 , hl2 = clean_bed_hlcalc(Kv, Ki, media_e, media_l, d, v_bed)
        hl3 = total_k*((v_pipe**2)/(2*g))
        h_total += (hl1+hl2+hl3)
        return h_total

    pump = st.selectbox('Pump selection',options=['30Z','55R','70R','70Z'])
    df_pump = pd.read_csv('Pumps.csv', header=[0,1])

    total_k = fitting_kvalue

    pumpx = df_pump[pump]['Flow']
    pumpy = df_pump[pump]['Head']
    sys_x = np.arange(0.1,30,0.5)
    sys_y = system_pressure(np.arange(0.1,30,0.5),total_k,valve_Kv,use_valve, use_darcy) +elev_hl

    pump_curve = LineString(np.column_stack((pumpx, pumpy)))
    system_curve = LineString(np.column_stack((sys_x, sys_y)))
    intersection = pump_curve.intersection(system_curve)

    #%% calculation for bed expansion


    plt.style.use('seaborn-dark')
    fig , ax = plt.subplots()
    ax.plot(sys_x, sys_y, label='System curve')
    ax.plot(pumpx,pumpy,label='Pump curve')
    ax.scatter(*intersection.xy,label='Duty point',color='orange',zorder=11)
    ax.set_xlim(left=0)
    ax.set_ylim(0,25)
    for expansion in [0, 20, 40, 60, 80, 100]:
        re_depth = media_l + (expansion / 100) * media_l
        re_e = 1 - ((media_l / re_depth) * (1 - media_e))

        res = back_wash_flow(Kv, Ki, d, re_e, media_rho)

        ax.axvline(res[0], linestyle='--',alpha=0.3, color='red')
        ax.text(res[0],21,str(expansion)+ '%' ,rotation=90)
    ax.set_xlabel('Flowrate $L/min$')
    ax.set_ylabel('Head $m$')
    ax.legend()
    ax.grid()
    st.pyplot(fig)
    col1, col2 = st.columns(2)
    col1.metric('Duty point Flow', str(np.around(intersection.x,4)) + ' L/min')
    col2.metric('Duty point Head', str(np.around(intersection.y,4)) + ' L/min')













