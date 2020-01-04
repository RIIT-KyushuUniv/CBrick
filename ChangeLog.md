# CBrick

## TODO

### diff3d
- MPIでの残差の値が異なる
- cell version



## REVISION HISTORY

---
- 2020-01-04  Version 1.4.3
  - intel_F_TCS >> ITO_TCS


---
- 2019-12-25 Version 1.4.2
  - add include <typeinfo> in CB_Comm_inline.h 


---
- 2019-03-06 Version 1.4.1
  - numPorc=1 の場合にも同じ記述ができるように、処理を追加
  
---
- 2019-03-06 Version 1.4.0
  - templated class
  - add int, unsigned, long long

---
- 2018-11-30 Version 1.3.2
  - node,  YZ方向チェック済
  - cellはコンパイルを通しただけで正しくない

---
- 2018-11-30 Version 1.3.1
  - 袖通信が怪しいので、テストを書く  > commtest
  - node, X方向のみ確認 >> 2, 3分割, F/Cindex
  - Y, Z方向の修正は、commtestの方向とpack/unpackを見ること


---
- 2018-11-10 Version 1.3.0
  - enumerate()でcellとnodeを共通化


---
- 2018-11-10 Version 1.2.2
  - CB_SubDomain.cpp (203):  bug fix if-statement
  - example のBrickCommクラス対応


---
- 2018-11-09 Version 1.2.1
  - getNumCandidates4IJ()
  - registerCandidates4IJ_Cell()
  - registerCandidates4IJ_Node()


---
- 2018-11-02 Version 1.2.0
  - 通信部分をComm classへ分離、IJK > KJIへ対応のため


---
- 2018-11-02 Version 1.1.1
  - beautify
  - rename functions >> packing


---
- 2018-10-02 Version 1.1.0
  - Bug fix: CB_CommS/V.cpp, MPI_Status stat[24] >> 26


---
- 2018-09-25 Version 1.0.6
  - checkOpenMPマクロを明示的に変更


---
- 2018-09-25 Version 1.0.5
  - Intel_F_TCS option


---
- 2018-09-24 Version 1.0.4
  - CMP0012
  - CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS TRUE


---
- 2018-09-11 Version 1.0.3
  - void getSubDomainHead(const int m, int m_sz[3])


---
- 2018-09-09 Version 1.0.2
  - void getSubDomainSize(const int m, int m_sz[3])


---
- 2018-03-18 Version 1.0.1
  - bool G2L_index(const int Gi, int* Li, const int c)


---
- 2018-03-03 Version 1.0.0
  1. 斜め方向通信の実装
    - エッジ、コーナーの斜め方向通信を実装。
    - 斜め方向通信を有効にする場合、コンパイルオプションに`-D_DIAGONAL_COMM`を追加
    - cmake実行時のオプション

  2. CB_Define.h
    - NOFACE、enum DIRectionの修正斜め方向フラグを追加（エッジ12、コーナー8）

  3. CB_SubDomain.h/cpp
    - (1) 隣接数配列サイズの変更
      - [6]で宣言されていた配列を[NOFACE]での宣言に変更（int cm[6] -> int cm[NOFACE]等
    ループもNOFACEで回すように修正

    - (2) 斜め方向袖通信バッファの宣言
      - エッジ、コーナー用の袖通信バッファの宣言を追加
      ~~~
      REAL_TYPE* f_es;   // edge send
      REAL_TYPE* f_er;   // edge recv
      REAL_TYPE* f_cs;   // corner send
      REAL_TYPE* f_cr;   // corner recv
      ~~~

    - (3) 斜め方向の通信テーブル作成（createRankTable関数）
    sd[m].cmに斜め方向（エッジ、コーナー）ランクを追加。

    - (4) 斜め方向袖通信関数の宣言追加
      - パックとisend/irecv
      ~~~
      bool pack_SE();  スカラー、エッジ、セル
      bool pack_SE(); スカラー、エッジ、ノード
      bool pack_SC();  スカラー、コーナー、セル
      bool pack_SC(); スカラー、コーナー、ノード
      bool pack_VE();  ベクター、エッジ、セル
      bool pack_VE(); ベクター、エッジ、ノード
      bool pack_VC();  ベクター、コーナー、セル
      bool pack_VC(); ベクター、コーナー、ノード
      ~~~

      - 展開
      ~~~
      void unpack_SE();  スカラー、エッジ、セル
      void unpack_SE(); スカラー、エッジ、ノード
      void unpack_SC();  スカラー、コーナー、セル
      void unpack_SC(); スカラー、コーナー、ノード
      void unpack_VE();  ベクター、エッジ、セル
      void unpack_VE(); ベクター、エッジ、ノード
      void unpack_VC();  ベクター、コーナー、セル
      void unpack_VC(); ベクター、コーナー、ノード
      ~~~

  4. CB_CommS.cpp
    - (1) initComm関数
      - 通常（面）袖通信のバッファサイズ変更
      - 通信方向の横方向の袖を通信しないよう調整
      ~~~
      f_sz[0] = size[1] * size[2] * gc * num_compo;
      f_sz[1] = size[0] * size[2] * gc * num_compo;
      f_sz[2] = size[0] * size[1] * gc * num_compo;
      ~~~

      - 斜め方向通信用送受信バッファの確保（initComm関数）

    - (2) Comm_S_nonblocking関数
      - req配列初期化ループ回転数の変更
        - [12]から[NOFACE*2]に変更

      - 通常（面）袖通信の通信サイズ変更
        - 通信方向の横方向の袖を通信しないよう調整
      ~~~
      msz[0] = size[1] * size[2] * gc_comm;
      msz[1] = size[0] * size[2] * gc_comm;
      msz[2] = size[0] * size[1] * gc_comm;
      ~~~

      - 斜め方向袖通信関数のコール
        - pack_SE()、pack_SEnode()、pack_SC()、pack_SCnode()をコール

    - (3) Comm_S_wait_nonblocking関数
      - wait用のMPI_Status配列のサイズ変更
      - [4]から[24]に変更（エッジの送受信時に12*2=24が最大数として必要なため）
      - 斜め方向袖通信のwaitと展開処理関数のコール追加
      - `unpack_SE(), unpack_SEnode(), unpack_SC(), unpack_SCnode()`をコール
      - (1)以外はCB_CommV.cppも同様の修正

  5. CB_Pack.h
    - (1) 面データのインデクス計算マクロの修正
      - 通信方向に横方向の袖を通信しないため、計算式を修正（_IDX_SI等）

  6. CB_PackingScalarCell.cpp
    - (1) 通常（面）袖通信の横方向袖の通信削除
      - 通信方向の横方向の袖を通信しないよう調整
      ~~~
      for( int k=0; k<kmax; k++ ){
      for( int j=0; j<jmax; j++ ){
      for( int i=0; i<vc_comm; i++ ){
      ~~~

    - (2) 斜め方向パック、通信、展開処理関数の実体追加
      - 上記3.(4)の処理関数を追加

      - 以下のソースファイルも同様の修正
      ~~~
      CB_PackingScalarNode.cpp
      CB_PackingVectorCell.cpp
      CB_PackingVectorNode.cpp
      ~~~

  7. example/diff3d/main.cpp
    - MPI_Request配列のサイズ変更
    - 配列サイズを[12]から[NOFACE*2]に変更
    - `MPI_Request req[NOFACE*2]`

  8. 袖通信テスト用のサンプルプログラムを追加
    - 以下に追加
    ~~~
    example/commtest
    example/CMakeFiles.txtも修正
    ~~~
   - 本サンプルプログラムを用いて、スカラー、ベクター、ノード、セルの斜め方向袖通信が正しく行われていることを確認。


---
- 2018-02-08 Version 0.9.11
  - bug fix:CB_Define.hとCB_Pack.hのマクロの修正 _VC >> (_VC)
    1. CB_Define.h
     - 3次元->1次元インデクス変換マクロの引数が計算式で渡された際に、正しく計算されない問題を修正。_IDX_S3D等
    2. CB_Pack.h
      - バッファへのインデクス変換マクロの引数が計算式で渡された際に、正しく計算されない問題を修正。_IDX_SI等

---
- 2018-01-31 Version 0.9.10
  - diff3dのfloatとdoubleのモジュールビルドは、現在のCMakeLists.txtでは、floatのC++オブジェクトにdoubleのFortranオブジェクトがリンクされてしまう。floatのみをビルドする。doubleはfloatをコメントアウトすればよいはず。

---
- 2018-01-30 Version 0.9.9
  - Split packing method into Cell and Node
  - bug fix : pack/unpack for node
  - modify SubDomain::setDivision(), div2.cpp, div4.cpp
  - bug fix : sortVolume
  - bug fix : counter in registerCandidates*()


---
- 2018-01-29 Version 0.9.8
  - expire serial build


---
- 2017-09-21 Version 0.9.7
  - Change notation of
    - X_minus >> I_minus, also others
    - `f_xms >> f_ims`, `f_xps >> f_ips`, `f_yms >> f_jms`,...
    - SX >> SI, SY >> SJ, SZ >> SK


---
- 2017-09-20 Version 0.9.6
  - bug fix : `registerCandidates4JK_Cell() / Node()` >>  tbl[odr].dsz[0]=G_size[0]


---
- 2017-09-04 Version 0.9.5
  - bug fix : `registerCandidates() >> registerCandidates_Cell() / Node()`
  - `registerCandidates4JK_Cell() / Node()`


---
- 2017-09-04 Version 0.9.4
  - virtual instruction modifier for destructor
  - JK division

---
- 2017-08-31 Version 0.9.3
  - use `CMAKE_BUILD_TYPE` for Debug > > define -DNDEBUG

---
- 2017-08-30 Version 0.9.2
  - change implementation of createRankTable() for OpenMP
  - add omp single

---
- 2017-08-30 Version 0.9.1
  - change instance method of cntl_tbl class
  - add buf_flag for delete pointer

---
- 2017-08-24 Version 0.9.0
  - Vector pack/unpack
  - remove blocking communication


---
- 2017-08-08 Version 0.8.0
  - build both float and double libraries


---
- 2017-08-03 Version 0.7.4
  - add stub for MPI_Comm and MPI_Request


---
- 2017-08-03 Version 0.7.3
  - add #ifndef DISABLE_MPI for serial build of applications >> expired later


---
- 2017-08-03 Version 0.7.2
  - found bug in blocking communication method, open issue.


---
- 2017-08-03 Version 0.7.1
  - expire use_NB


---
- 2017-08-03 Version 0.7.0
  - modify G2L_index()


---
- 2017-08-03 Version 0.6.9
  - modify getHeadIndex()


---
- 2017-08-02 Version 0.6.8
  - write size, head info into div_process.txt


---
- 2017-08-02 Version 0.6.7
  - bug fix of moemory release in findOptimalDivision() and findParameter()


---
- 2017-07-28 Version 0.6.6
  - change scope of Comm_V_blocking() to public


---
- 2017-07-28 Version 0.6.5
  - setSubDomain()


---
- 2017-07-28 Version 0.6.4
  - useNonblocking()


---
- 2017-07-28 Version 0.6.3
  - CB_Version.h


---
- 2017-07-28 Version 0.6.2
  - bug fix : if ( !f ) delete[] f; => if (f) ...


---
- 2017-07-28 Version 0.6.1
  - fix install doc file


---
- 2017-07-27 Version 0.6.0
  - install header files


---
- 2017-07-27 Version 0.5.0
  - convert to Fortran index >> f_index


---
- 2017-07-27 Version 0.4.0
  - example/diff3d
  - add Fortran env in case of turning on example
  - introduce auto_div, and merge findParameter() into findOptimalDivision()


---
- 2017-07-26 Version 0.3.0
  - pack/unpack
  - CB_CommV.cpp
  - sortLenX(), sortCube()
  - findParameter()
  - example/div3.cpp, div4.cpp


---
- 2017-07-23 Version 0.2.1
  - CB_Comm.cpp


---
- 2017-07-22 Version 0.2.0
  - createRankTable()

---
- 2017-07-20 Version 0.1.0
  - Cmake compilation
  - Node and Cell


---
- 2017-06-29 Version 0.0.1
  - initial material
  - serial compile by makefile
    - `CB_SubDomain.cpp` and `example/division/div.cpp`
