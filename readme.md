# README.MD for ANALYSIS OF HST OBSERVATION

## To-Do List
- ✅ 木星が視野に入っているSTIS観測のObs. IDをリストアップする。
- プロポーザルを読み、おおよその視野を調べる。
- Campaign中のCMLを調べて、IDの絞り込みをする。
    - ✅ UTCと赤緯赤経はリストに含まれる。
    - System IIIはSpiceypyで計算できる？
        1. ✅ SpiceypyとAstropyのインストール。
        1. Spiceypyの使い方を勉強する。
        1. JUICE用のモジュールがある。打ち上げ前の時刻にも対応してる？木星から見た衛星の位置も計算可能。
            - "For example, you can get the Io position compared to Jupiter in the 'J2000' frame."
            - https://irfpy.irf.se/projects/pep/issues/issue_juice_spice_readme.html
        1. 各衛星のSys. III経度を求める。