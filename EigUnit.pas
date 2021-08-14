unit EigUnit;

interface

uses
  DAVLib.Types;

function eig(A, B: TDComplexDynMatrix): TDComplexDynArray; overload;

function deig(A, B: TDoubleDynMatrix): TDComplexDynArray;

implementation

uses
  System.SysUtils, System.Math,  DAVLib.Lapack, DAVLib.Blas,  DAVLib.ComplexType, DAVLib.SimpleTypesHelper,
  DAVLib.MatrixTypesHelper, DAVLib.VectorTypesHelper, DAVLib.ComplexMatrixTypesHelper, DAVLib.ComplexVectorTypesHelper;

function eig(A, B: TDComplexDynMatrix): TDComplexDynArray;
var
  BALANC, JOBVL, JOBVR, SENSE: Char;
  N: Integer;
  LDA: Integer;
  LDB: Integer;
  ALPHA, BETA, VL, VR, WORK: TDComplexDynArray;
  LSCALE, RSCALE, RCONDE, RCONDV, RWORK: TDoubleDynArray;
  ILO, IHI, INFO: Integer;
  ABNRM, BBNRM: Double;
  IWORK: TIntegerDynArray;
  BWORK: TBooleanDynArray;
  LDVL: Integer;
  LDVR: Integer;
  LWORK, LRWORK: Integer;
  I, J, K: Integer;
begin
  BALANC := 'N';
  JOBVL := 'N';
  JOBVR := 'N';
  SENSE := 'E';
  N := A.RowCount; // col for fortran column-major order
  LDA := N;
  LDB := N;
  LDVL := 1;
  LDVR := 1;
  LWORK := 2 * N * N + 2 * N;
  LRWORK := 6 * N;

  ALPHA := TDComplexDynArray.Create(N);
  BETA := TDComplexDynArray.Create(N);

  VL := TDComplexDynArray.Create(LDVL * N);
  VR := TDComplexDynArray.Create(LDVR * N);

  LSCALE := TDoubleDynArray.Create(N);
  RSCALE := TDoubleDynArray.Create(N);
  RCONDE := TDoubleDynArray.Create(N);
  RCONDV := TDoubleDynArray.Create(N);

  WORK := TDComplexDynArray.Create(LWORK);
  RWORK := TDoubleDynArray.Create(LRWORK);
  IWORK := TIntegerDynArray.Create(N + 2);
  SetLength(BWORK, N);

  ILO := 1;
  IHI := N;
  INFO := 0;
  ABNRM := 0;
  BBNRM := 0;

  try
    // 'Transposed' for fortran column-major order
    zggevx_(@BALANC, @JOBVL, @JOBVR, @SENSE, //
      @N, A.Transposed.ToVector.PItem[0], @LDA, //
      B.Transposed.ToVector.PItem[0], @LDB, //
      ALPHA.PItem[0], BETA.PItem[0], //
      VL.PItem[0], @LDVL, VR.PItem[0], @LDVR,  //
      @ILO, @IHI, LSCALE.PItem[0], RSCALE.PItem[0], //
      @ABNRM, @BBNRM, RCONDE.PItem[0], RCONDV.PItem[0], //
      WORK.PItem[0], @LWORK, RWORK.PItem[0], IWORK.PItem[0], @BWORK[0], @INFO);
  finally
    if INFO <> 0 then
      raise Exception.Create('Lapack: Error ' + INFO.ToString);
  end;

  Result := TDComplexDynArray.Create(ALPHA);
  for I := 0 to N - 1 do
  begin
    if not BETA[I].IsZero then
      Result[I] := Result[I] / BETA[I];
  end;
end;

function deig(A, B: TDoubleDynMatrix): TDComplexDynArray;
var
  BALANC, JOBVL, JOBVR, SENSE: Char;
  N: Integer;
  LDA: Integer;
  LDB: Integer;
  ALPHAR, ALPHAI, BETA, VL, VR, WORK: TDoubleDynArray;
  LSCALE, RSCALE, RCONDE, RCONDV: TDoubleDynArray;
  ILO, IHI, INFO: Integer;
  ABNRM, BBNRM: Double;
  IWORK: TIntegerDynArray;
  BWORK: TBooleanDynArray;
  LDVL: Integer;
  LDVR: Integer;
  LWORK: Integer;
  I, J, K: Integer;
begin
  BALANC := 'N';
  JOBVL := 'N';
  JOBVR := 'N';
  SENSE := 'E';
  N := A.ColCount; // col for fortran column-major order
  LDA := N;
  LDB := N;
  LDVL := 1;
  LDVR := 1;
  LWORK := 2 * N * N + 2 * N;

  ALPHAR := TDoubleDynArray.Create(N);
  ALPHAI := TDoubleDynArray.Create(N);
  BETA := TDoubleDynArray.Create(N);

  VL := TDoubleDynArray.Create(LDVL * N);
  VR := TDoubleDynArray.Create(LDVR * N);

  LSCALE := TDoubleDynArray.Create(N);
  RSCALE := TDoubleDynArray.Create(N);
  RCONDE := TDoubleDynArray.Create(N);
  RCONDV := TDoubleDynArray.Create(N);

  WORK := TDoubleDynArray.Create(LWORK);
  IWORK := TIntegerDynArray.Create(N + 2);
  SetLength(BWORK, N);

  ILO := 1;
  IHI := N;
  INFO := 0;
  ABNRM := 0;
  BBNRM := 0;

  try
    // 'Transposed' for fortran column-major order
    dggevx_(@BALANC, @JOBVL, @JOBVR, @SENSE, //
      @N, A.Transposed.ToVector.PItem[0], @LDA,
      B.Transposed.ToVector.PItem[0], @LDB, //
      ALPHAR.PItem[0], ALPHAI.PItem[0], BETA.PItem[0], //
      VL.PItem[0], @LDVL, VR.PItem[0], @LDVR,  //
      @ILO, @IHI, LSCALE.PItem[0], RSCALE.PItem[0], //
      @ABNRM, @BBNRM, RCONDE.PItem[0], RCONDV.PItem[0], //
      WORK.PItem[0], @LWORK, IWORK.PItem[0], @BWORK[0], @INFO);
  finally
    if INFO <> 0 then
      raise Exception.Create('Lapack: Error ' + INFO.ToString);
  end;

  Result := TDComplexDynArray.Create(ALPHAR, ALPHAI);
  for I := 0 to N - 1 do
  begin
    if not BETA[I].IsZero then
      Result[I] := Result[I] / BETA[I];
  end;
end;

end.

