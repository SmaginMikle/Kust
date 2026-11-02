# Экспорт результатов Python-двойника в тот же JSON/CSV/HTML, что формирует C# (Reporting.cs).
import json, math, os, sys, pickle, datetime
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from solver import P, Terrain, Problem, pad_len

inp, mapf, outdir, tpl = sys.argv[1:5]
(results, infeasible) = pickle.load(open(os.path.join(outdir, 'results.pkl'), 'rb'))
terr = Terrain(mapf); prob = Problem(inp, terr)
W = prob.wells
names = {'1': 'Суходол', '2': 'Болото', '3': 'Водоём'}

def r1(x, d=1):
    return None if x is None or (isinstance(x, float) and not math.isfinite(x)) else round(float(x), d)

# лучший — допустимый оптимизационный вариант первого сценария
cands = [r for r in results if r['kind'] == 'opt' and not r.get('same_as') and r.get('scenario', '').startswith('Разворот')]
best = min(cands, key=lambda r: (r['violations'], r['total_cost']))

# участки (как Reporting.Regions)
n = len(W); lab = list(range(n))
def find(x):
    while lab[x] != x: lab[x] = lab[lab[x]]; x = lab[x]
    return x
for i in range(n):
    for k in range(i + 1, n):
        if math.hypot(W[i]['T1'][0] - W[k]['T1'][0], W[i]['T1'][1] - W[k]['T1'][1]) < 5000: lab[find(i)] = find(k)
groups = {}
for i in range(n): groups.setdefault(find(i), []).append(i)
groups = sorted(groups.values(), key=lambda g: np.mean([W[i]['T1'][0] for i in g]))
regions = []; m = P.MaxOffset + 300
for no, g in enumerate(groups, 1):
    e0 = min(min(W[i]['T1'][1], W[i]['T3'][1]) for i in g) - m; e1 = max(max(W[i]['T1'][1], W[i]['T3'][1]) for i in g) + m
    n0 = min(min(W[i]['T1'][0], W[i]['T3'][0]) for i in g) - m; n1 = max(max(W[i]['T1'][0], W[i]['T3'][0]) for i in g) + m
    e0, n0, e1, n1 = math.floor(e0 / 100) * 100, math.floor(n0 / 100) * 100, math.ceil(e1 / 100) * 100, math.ceil(n1 / 100) * 100
    c0 = max(0, int(math.floor((e0 - terr.x0) / terr.cs))); c1 = min(terr.nc - 1, int(math.floor((e1 - terr.x0) / terr.cs)))
    rr0 = max(0, terr.nr - 1 - int(math.floor((n1 - terr.y0) / terr.cs))); rr1 = min(terr.nr - 1, terr.nr - 1 - int(math.floor((n0 - terr.y0) / terr.cs)))
    grid = terr.g[rr0:rr1 + 1, c0:c1 + 1].copy(); grid[grid == terr.nodata] = 9
    regions.append(dict(name=f"Участок {no} ({len(g)} скв.)", bounds=[e0, n0, e1, n1], gridE0=terr.x0 + c0 * terr.cs,
                        gridN1=terr.y0 + (terr.nr - rr0) * terr.cs, cell=terr.cs, grid=grid.tolist()))

def variant_json(r):
    d = dict(name=r['name'], kind=r['kind'], scenario=r.get('scenario', 'Исходные данные'), K=r.get('K'), sameAs=r.get('same_as'),
             nPads=r['n_pads'], nWells=r['n_wells'], mdTotal=r1(r['md_total']), drillCost=r1(r['drill_cost']), padBuildCost=r1(r['pad_build_cost']),
             sandCost=r1(r['sand_cost']), sandVolume=r1(r['sand_volume']), totalCost=r1(r['total_cost']), violations=r['violations'],
             estimate=r1(r.get('opt_estimate')))
    if r.get('same_as'): return d
    d['pads'] = [dict(name=p['name'], n=r1(p['N'], 2), e=r1(p['E'], 2), count=p['n'], angle=r1(p['angle']), length=r1(p['length']), width=r1(p['width']),
                      terrain={k: round(v, 3) for k, v in p['terrain'].items()}, sandVolume=r1(p['sand_volume'], 0), sandCost=r1(p['sand_cost']),
                      buildCost=r1(p['build_cost']), forbidden=p['forbidden'], wells=p['wells']) for p in r['pads']]
    ws = []
    slot = {}
    for w in r['wells']:
        slot[w['pad']] = slot.get(w['pad'], 0) + 1
        sv = w.get('survey')
        tr = []
        if sv is not None:
            for i, s in enumerate(sv):
                if i % 5 == 0 or i == len(sv) - 1 or abs(s[0] - w['mdT1']) < 1e-6:
                    tr.append([round(s[3], 1), round(s[4], 1), round(s[5], 1), round(s[0], 1)])
        warn = [f"сближение {w['min_sep']:.1f} м"] if w['min_sep'] is not None and w['min_sep'] < 5 else []
        ws.append(dict(well=w['well'], pad=w['pad'], slot=slot[w['pad']], whN=r1(w['whN'], 2), whE=r1(w['whE'], 2), offset=r1(w['offset']),
                       azTurn=r1(w['azturn']), hz=r1(w['hz']), mdT1=r1(w['mdT1']), mdT3=r1(w['mdT3']), tanLen=r1(w.get('tangent')),
                       tanInc=r1(w.get('inc_tangent'), 2), tanAz=r1(w.get('az_tangent'), 2), cost=r1(w['cost']), costNNS=r1(w.get('cost_nns')),
                       costGS=r1(w.get('cost_gs')), hook=r1(w['hook']), torque=r1(w['torque']), maxDls=r1(w.get('maxdls'), 3),
                       minSep=r1(w['min_sep']), viol=w['viol'], warn=warn, tr=tr))
        w['warn'] = warn; w['slot'] = slot[w['pad']]
    d['wells'] = ws
    return d

R = dict(meta=dict(input='input_mass.json', map=os.path.basename(mapf), generated=datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
                   engine='Python-двойник WellPadOptimizer', best=best['name'], infeasible=infeasible,
                   params=dict(MaxOffset=P.MaxOffset, MaxMD=P.MaxMD, Dls1=P.Dls1, Dls2=P.Dls2, VerticalLength=P.VertLen, UseAzimuthTurn=True,
                               MaxAzimuthTurn=P.MaxAzTurn, MaxHookLoad=P.MaxHookLoad, MaxTorque=P.MaxTorque, CostNNS=P.CostNNS, CostGS=P.CostGS,
                               SandCostPer1000m3=P.SandCostPer1000, MinWellsPerPad=P.MinWells, StartMaxWellsPerPad=P.StartMaxWells, EndMaxWellsPerPad=P.EndMaxWells)),
         forbiddenCode=9, terrainClasses=[dict(code=c, name=names[str(c)], fill=h) for c, h in P.FillHeight.items()],
         regions=regions, wells=[dict(name=w['name'], uid=w['uid'], t1=[round(x, 2) for x in w['T1']], t3=[round(x, 2) for x in w['T3']]) for w in W],
         variants=[variant_json(r) for r in results])
js = json.dumps(R, ensure_ascii=False, separators=(',', ':'))
open(os.path.join(outdir, 'results.json'), 'w', encoding='utf-8').write(js)
body = open(tpl, encoding='utf-8').read().replace('/*__REPORT_DATA__*/null', js.replace('</', '<\\/'))
open(os.path.join(outdir, 'report_fragment.html'), 'w', encoding='utf-8').write(body)
open(os.path.join(outdir, 'report.html'), 'w', encoding='utf-8').write(
    '<!doctype html>\n<html lang="ru">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
    '<style>body{margin:0}[hidden]{display:none!important}</style>\n</head>\n<body>\n' + body + '\n</body>\n</html>\n')

# ---- CSV (как C#)
def F(v, d=1):
    if v is None or (isinstance(v, float) and not math.isfinite(v)): return ''
    return f"{v:.{d}f}".replace('.', ',')
bom = '﻿'
L = ['Вариант;Сценарий;Макс. скв. на кусте;Совпадает с;Кустов;Скважин;Сумма MD, м;Бурение, тыс.руб;Строительство кустов, тыс.руб;Отсыпка, м3;Отсыпка, тыс.руб;Итого, тыс.руб;Нарушений']
for r in results:
    L.append(';'.join(map(str, [r['name'], r.get('scenario', 'Исходные данные'), r.get('K') or '', r.get('same_as') or '', r['n_pads'], r['n_wells'],
                                 F(r['md_total']), F(r['drill_cost']), F(r['pad_build_cost']), F(r['sand_volume'], 0), F(r['sand_cost']), F(r['total_cost']), r['violations']])))
for x in infeasible:
    L.append(';'.join(map(str, [x['name'], x['scenario'], x['K'], '', '', '', '', '', '', '', '', 'нет решения: ' + x['reason'], ''])))
open(os.path.join(outdir, 'summary.csv'), 'w', encoding='utf-8').write(bom + '\n'.join(L) + '\n')
L = ['Вариант;Куст;N;E;Скважин;Угол ряда, °;Длина, м;Ширина, м;Суходол, доля;Болото, доля;Водоём, доля;Отсыпка, м3;Отсыпка, тыс.руб;Строительство, тыс.руб;Запретная зона;Скважины']
for r in results:
    if r.get('same_as'): continue
    for p in r['pads']:
        L.append(';'.join(map(str, [r['name'], p['name'], F(p['N'], 2), F(p['E'], 2), p['n'], F(p['angle'], 0), F(p['length']), F(p['width']),
                                     F(p['terrain']['1'], 3), F(p['terrain']['2'], 3), F(p['terrain']['3'], 3), F(p['sand_volume'], 0), F(p['sand_cost']), F(p['build_cost']),
                                     'да' if p['forbidden'] else 'нет', ' '.join(p['wells'])])))
open(os.path.join(outdir, 'pads.csv'), 'w', encoding='utf-8').write(bom + '\n'.join(L) + '\n')
L = ['Вариант;Скважина;UID;Куст;Слот;Устье N;Устье E;T1 N;T1 E;T1 TVD;T3 N;T3 E;T3 TVD;Отход до T1, м;Разворот, °;MD до T1, м;MD до T3, м;ГС, м;Тангенс, м;Зенит тангенса, °;Азимут тангенса, °;Макс. ПИ, °/10м;Нагрузка на крюке, тс;Момент, кН·м;Мин. сближение, м;ННС, тыс.руб;ГС, тыс.руб;Бурение, тыс.руб;Нарушения;Предупреждения']
for r in results:
    if r.get('same_as'): continue
    for w in r['wells']:
        L.append(';'.join(map(str, [r['name'], w['well'], w['uid'], w['pad'], w['slot'], F(w['whN'], 2), F(w['whE'], 2), F(w['T1'][0], 2), F(w['T1'][1], 2), F(w['T1'][2], 2),
                                     F(w['T3'][0], 2), F(w['T3'][1], 2), F(w['T3'][2], 2), F(w['offset']), F(w['azturn']), F(w['mdT1']), F(w['mdT3']), F(w['hz']),
                                     F(w.get('tangent')), F(w.get('inc_tangent'), 2), F(w.get('az_tangent'), 2), F(w.get('maxdls'), 2), F(w['hook']), F(w['torque']),
                                     F(w['min_sep']), F(w.get('cost_nns')), F(w.get('cost_gs')), F(w['cost']), ', '.join(w['viol']), ', '.join(w['warn'])])))
open(os.path.join(outdir, 'wells.csv'), 'w', encoding='utf-8').write(bom + '\n'.join(L) + '\n')
os.makedirs(os.path.join(outdir, 'trajectories'), exist_ok=True)
for vi, r in enumerate(results, 1):
    if r.get('same_as'): continue
    L = ['Скважина;Куст;MD, м;Зенит, °;Азимут, °;N;E;TVD, м;ПИ, °/10м']
    for w in r['wells']:
        if w.get('survey') is None: continue
        for s in w['survey']:
            L.append(';'.join([w['well'], w['pad'], F(s[0], 2), F(s[1], 3), F(s[2], 3), F(s[3], 2), F(s[4], 2), F(s[5], 2), F(s[6], 3)]))
    open(os.path.join(outdir, 'trajectories', f'variant_{vi:02d}.csv'), 'w', encoding='utf-8').write(bom + '\n'.join(L) + '\n')
print('best', best['name'], round(best['total_cost'] / 1000, 1), 'json', len(js) // 1024, 'KB')
