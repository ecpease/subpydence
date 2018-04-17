import numpy as np



enums = [1,1,1]

dhcEnum, vkEnum, sfEnum = enums

def getCompParams(vkEnum, sfEnum, depth):  #

    # Enum   0 = opt, 1 = avg, 2 = pess

    ########vk
    print(depth)
    depth = max(depth, 1.0)  # no negative depths

    vk_min = 1.02E-5 * m.exp((-2.55E-3) * depth)
    vk_max = 0.0342 * (depth ** (-1.321))
    vk_avg = (vk_min + vk_max) / 2

    if vkEnum == 0:
        vk = vk_min  # smaller vk, less subsidence
    elif vkEnum == 1:
        vk = vk_avg
    elif vkEnum == 2:
        vk = vk_max
    else:
        raise ValueError('Enum issues')
    #######sfv, sfe
    por = 1.4485 * (depth ** -0.233)
    comp = 0.000003 * (depth ** -0.703)
    scv = 1000 * 9.81 * comp / 3.2808
    sce = scv * 0.01

    sfv_avg = scv + 1000 * 9.81 * por * 4.4E-10 / 3.2808
    sfv_min = sfv_avg * 0.3
    sfv_max = sfv_avg * 3.3

    sfe_avg = sce + 1000 * 9.81 * por * 4.4E-10 / 3.2808
    sfe_min = sfe_avg * 0.3
    sfe_max = sfe_avg * 3.3

    if sfEnum == 0:
        sfv = sfv_min  # small compress, less subsidence
        sfe = sfe_min
    elif sfEnum == 1:
        sfv = sfv_avg
        sfe = sfe_avg
    elif sfEnum == 2:
        sfv = sfv_max
        sfe = sfe_max
    else:
        raise ValueError('Enum issues')

    return vk, sfe, sfv


def GetComposite(thicks):
    sumBsq = 0
    sumB = 0
    ndb = 0
#    print(thicks)
    for db in thicks:
        ndb+=1
        sumBsq += db ** 2
        sumB += db
    bequiv = np.sqrt(sumBsq / ndb)
    nequiv = sumB / bequiv
    dz = bequiv
    rnb = nequiv
    return (dz,rnb)






def get_sub_params(dis,ssHead,enums,nlay):
    for mlay in model_layers:  # iterate by "layer"
        print(mlay)
        dzAry = []
        rnbAry = []
        # gamxAry = []
        # gamyAry = []
        # apiSub = {}
        # apiMlay = {}
        # apiRC = {}
        rAry = []
        cAry = []
        apiAry = []
        nmzAry = []
        vkAry = []
        if (l > 0):  # skip layer 1
            for api in apis:
                #        for api in apis[0:15]:
                df_clay_thick_af = df_clay_thick.loc[
                    (df_clay_thick['API_10'] == api) & (df_clay_thick['model_layer'] == mlay)]
                if (not df_clay_thick_af.empty):
                    thicks = df_clay_thick_af['Clay_length']
                    (dz, rnb) = GetComposite(thicks)
                    gamx = df_clay_thick_af['xGAM'].iloc[0]
                    gamy = df_clay_thick_af['yGAM'].iloc[0]
                    r, c = GetRCfromXY(x=gamx, y=gamy, x0=dis.sr.xul, y0=dis.sr.yul, ang=dis.sr.rotation, delR=delR,
                                       delC=delC)  # one based r,c
                    offset = 0
                    while (r == 0 or c == 0):
                        offset += 1000
                        ang = dis.sr.rotation
                        xnew = gamx + offset * m.cos(ang * np.pi / 180)
                        ynew = gamy + offset * m.sin(ang * np.pi / 180)
                        r, c = GetRCfromXY(x=xnew, y=ynew, x0=dis.sr.xul, y0=dis.sr.yul, ang=dis.sr.rotation, delR=delR,
                                           delC=delC)
                        #                        print "moving r and c to get a better depth using fixed angle of 36.7", api, r, c, xnew, ynew

                    if (l == 1):  # pumping layer <---------- subject to change
                        thk1 = thickness[l][r - 1, c - 1]
                        thk2 = thickness[l + 1][r - 1, c - 1]
                        rnb *= thk1 / (thk1 + thk2)
                    elif (l == 2):  # layer below pumping layer
                        thk1 = thickness[l - 1][r - 1, c - 1]
                        thk2 = thickness[l][r - 1, c - 1]
                        rnb *= thk2 / (thk1 + thk2)

                    # thk1 = thickness[l - 1][r - 1, c - 1]
                    # thk2 = thickness[l][r - 1, c - 1]
                    # rnb *= thk2 / (thk1 + thk2)

                    dzAry.append(dz)
                    rnbAry.append(rnb)
                    # apiSub[api+" "+form] = (dz,rnb)
                    apiMlay[api + '_' + str(mlay)] = (dz, rnb)
                    rAry.append(r)
                    cAry.append(c)

                    apiAry.append(api)
                    gamxAry.append(gamx)
                    gamyAry.append(gamy)
                    ndb += 1
                    nmz += 1
                    nmzapi[api] = nmz
                    #                  if(nmz == 123):
                    #                       print "r,c,123",r,c
                    vk, sfe, sfv = getCompParams(vkEnum, sfEnum, depth[l, r - 1, c - 1])  # zero based depth
                    dp.append([vk, sfe, sfv])
                    nmzAry.append(nmz)
                    vkAry.append(vk)

                    ldn.append(l)

                    #        if(l == 0):
                    #            df = pd.DataFrame({'api':apiAry,'gamx':gamxAry, 'gamy':gamyAry, "nz":nmzAry, "vk":vkAry})
                    #            print "debug_enum", dhcEnum,vkEnum,sfEnum
                    #
                    #            df.to_csv(os.path.join('sub_api_diag_{}_{}_{}.csv'.format(dhcEnum,vkEnum,sfEnum)))
                    #            exit()

                    #        outDF = pd.DataFrame({'API':apiAry,'GAMx':gamxAry,'GAMy':gamyAry,'dz':dzAry,'rnb':rnbAry,'r':rAry,'c':cAry})
                    #        outDF.to_csv(os.path.join('smaller_model','sub_figs','dz_rnb_sub_'+form+'.csv'),index=False)

                    # we don't know for each formation how many ndbs we are adding, so resize array on the fly
        if (l > 0):  # late game skip Burkeville
            #            if(form == forms[0]):
            if (l == 1):
                rnb3d = np.zeros([ndb, nrow, ncol])
                dz3d = np.zeros([ndb, nrow, ncol])
                nz3d = np.zeros([ndb, nrow, ncol])
                dhc3d = np.zeros([ndb, nrow, ncol])
                dstart3d = np.zeros([ndb, nrow, ncol])

            else:
                rnb3d.resize(ndb, nrow, ncol)
                dz3d.resize(ndb, nrow, ncol)
                nz3d.resize(ndb, nrow, ncol)
                dhc3d.resize(ndb, nrow, ncol)
                dstart3d.resize(ndb, nrow, ncol)

            for r in range(nrow):
                for c in range(ncol):
                    x, y = GetXYfromRC(r + 1, c + 1, x0=dis.sr.xul, y0=dis.sr.yul, ang=dis.sr.rotation, delR=delR,
                                       delC=delC)
                    api = GetNearest(x, y, apiAry, gamxAry, gamyAry)
                    apiRC[(r, c)] = api

            for api in apiAry:
                for r in range(nrow):
                    for c in range(ncol):
                        apiClose = apiRC[(r, c)]
                        if (api == apiClose):
                            # dz,rnb = apiSub[api+" "+form]
                            dz, rnb = apiMlay[api + '_' + str(mlay)]
                            rnb3d[dbNum, r, c] = rnb
                            dz3d[dbNum, r, c] = dz
                            nz3d[dbNum, r, c] = nmzapi[api]
                            depthc = depth[ldn[dbNum], r, c]
                            ssHeadc = ssHead[ldn[dbNum], r, c]
                            dhc3d[dbNum, r, c] = getDHC(dhcEnum, depthc, ssHeadc)
                            #                            print l, r, c, dbNum, dhc3d[dbNum,r,c]
                            dstart3d[dbNum, r, c] = ssHead[ldn[dbNum], r, c]

                dbNum += 1
                #       if(l > 0):
                #           print "dhc", l, dhc3d[0,0,10]
        l += 1  # increment layer

    return (ndb, nmz, ldn, rnb3d, dp, dstart3d, dhc3d, dz3d, nz3d)


def compaction(ln, hc, sfe, sfv, com = 0):
    '''

    :param ln:  is a one-dimensional array specifying the model layer assignments for each system of no-delay interbeds.
    :param hc: HC is an array specifying the preconsolidation head or preconsolidation stress in terms of head in the
            aquifer for systems of no-delay interbeds. For any model cells in which specified HC is greater than the
            corresponding value of starting head, the value of HC will be set to that of starting head.
    :param sfe: Sfe is an array specifying the dimensionless elastic skeletal storage coefficient for systems of
            no-delay interbeds.
    :param sfv: Sfv is an array specifying the dimensionless inelastic skeletal storage coefficient for systems of
            no-delay interbeds.
    :param com: COM is an array specifying the starting compaction in each system of no-delay interbeds. Compaction
            values computed by the package are added to values in this array so that printed or stored values of
            compaction and land subsidence may include previous components. Values in this array do not affect
            calculations of storage changes or resulting compaction. For simulations in which output values are to
            reflect compaction and subsidence since the start of the simulation, enter zero values for all elements of
            this array.
    :return:
    '''
    
    
    




