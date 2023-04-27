# README.MD for ANALYSIS OF HST OBSERVATION

## To-Do List
- ✅ 木星が視野に入っているSTIS観測のObs. IDをリストアップする。
- プロポーザルを読み、おおよその視野を調べる。
    - ✅ ID 13035 のプロポーザルを入手。
        - https://www.stsci.edu/hst/phase2-public/13035.pdf
    - ✅ CML rangeの記載あり。
    - "NOT OLG OF JUPITER FROM EARTH BETWEEN 15 345"の意味は?
- Campaign中のCMLを調べて、IDの絞り込みをする。
    - ✅ UTCと赤緯赤経はリストに含まれる。
    - 🚫 System IIIはSpiceypyで計算できる？
    - ✅ JPLのHorizons System (https://ssd.jpl.nasa.gov/horizons/app.html#/) で調べることができる。
        - (5) "Table Settings" で (14) "Observer sub-lon & sub-lat" ($\phi_m$) を選択すると、"sub-lon"はSystem III経度で表示される。
        - 衛星の場合、(26) "Observer-Primary-Target angle" ($\varphi_{obs}$)を追加する。
        - (Sub-observer Sys. III からの離角) $ = 180^\circ$ $- \varphi_{obs} - (360^\circ - \phi_m)$ でOK。
        - 光速のラグ(だいたい40分 参考 http://photon.sci-museum.kita.osaka.jp/question/text/distance.html)は考慮済み。
    - ✅ SpiceypyとAstropyのインストール。
        1. Spiceypyの使い方を勉強する。
        1. JUICE用のモジュールがある。打ち上げ前の時刻にも対応してる？木星から見た衛星の位置も計算可能。
            - "For example, you can get the Io position compared to Jupiter in the 'J2000' frame."
            - https://irfpy.irf.se/projects/pep/issues/issue_juice_spice_readme.html
            - J2000座標は出るけど、Sys IIIは分からない。
        1. 各衛星のSys. III経度を求める。
- ✅ サンプル画像の表示　`sample_viewer.ipynb`
    - Nicholsさんから頂いたサンプル1枚。データIDは `OEOW08FGQ` 。
    - Fitsファイル単体の表示はOK。
    - 極の位置を特定できていない。なので緯度経度の表示もまだ。
    - System III経度は、UTから計算?
- 磁場モデルからFootprintのL値を調べる。
    - ❌ ルジャンドル陪関数を計算する `pyshtools` がARM64に対応していない。`conda install -c conda-forge pyshtools`はApple Siliconでは不可。
    - Fortran95用の `shtools` はApple Siliconでも動くらしい。でもmoduleの使い方分からない。
    - ✅ Deimosで磁場モデルを扱えるように環境を整備した。
- ✅ 各IDの積分時間の取り扱い。
    - `TIME-TAG`の説明 (https://hst-docs.stsci.edu/stisdhb/chapter-2-stis-data-structure/2-3-stis-file-structures)
    - `astropy.table`を使って`oeow08fgq_tag.fits`の中身を見ると、トータル積分時間に比べて細かい時間刻みで記録が格納されている。
        - Time-tagデータは、トータルの露光時間の範囲で任意の積分時間を指定できる。Nicholsさんは30sを指定している。
        - refデータにtime-tagが格納されている。
    - `oeow08fgq_raw.fits`から各時刻のデータが取り出せる？
- 座標系の取り扱い。
    - 1次元のflattenedデータを緯度経度マップにマッピングするコードを頂いた。
        - ❌ Errorが出て今のところうまく動かない。
- Multiple Spots of EFP Aurora
    - DetectionはAug 15 1999 by STIS (南半球4例、北半球1例) (Bonfond+2017)
    - 2022-09-28T21:25:22 (DoY 271)に2点のスポットを確認 (北半球)。