%%% @doc
%%%
%%%
%%% @author Ahmed Omar <spawn.think@gmail.com>
-module(matrices).
%%% @end
-define(DEFAULT, mat).
-export([new/0,new/1,
         delete/0,delete/1,
         from_list/2,to_list/1,
         lookup/2,update/2,foldl/3]).
-export([get_size/1,min/1,max/1,min_uncov/3]).
-export([is_covered/2,covered_cols/2,
         zero_in_rows/4,star_in_col/3,
         star_in_row/3,prime_in_row/3]).
-export([cover/2,uncover/2,uncover_all/1,
         prime/2,erase_prime/2,erase_primes/1,
         convert_path/3,star_zeros/5,
         sub_min/2]).


new()->
    new(?DEFAULT).
new(M)->
    ets:new(M,[named_table]).
delete()->
    delete(?DEFAULT).
delete(M)->
    ets:delete(M).

from_list(M, List)->
    ets:insert(M, List).
to_list(M)->
    ets:tab2list(M).

lookup(M,Key)->
    case ets:lookup(M,Key) of
        [{Key, Val}]->
            {Key,Val};
        _ ->
            false
    end.
update(M,E)->
    true = ets:insert(M,E).

get_size(M)->
   proplists:get_value(size,ets:info(M)).

foldl(Fun, Acc, M)->
    ets:foldl(Fun, Acc,M).

star(M,Key) ->
    update(M, {Key,1}).
prime(M,Key)->
    update(M, {Key,2}).
erase_prime(M, Key)->
    update(M, {Key,0}).
cover(M,Key)->
    update(M, {Key,true}).
uncover(M,Key)->
    update(M, {Key,false}).
convert_path(M, Path,I)->
    convert_path_h(M, Path, 0, I).
convert_path_h(_M,_Path,I,I)->
    ok;
convert_path_h(M, Path, X,_I)->
    {_,R}=lookup(Path,{X,0}),
    {_,C}=lookup(Path,{X,1}),
    {_,MVal}=lookup(M,{R,C}),
    update(M,{{R,C},convert(MVal)}),
    convert_path_h(M,Path,X+1,_I).
convert(1)->
    0;
convert(_) ->
    1.
uncover_all(M)->
    foldl(fun({Key,_},_Acc)->
                  uncover(M, Key),
                  _Acc;
             (_,_Acc) ->
                  _Acc
          end,[], M).

erase_primes(M)->
    foldl(fun({Key, 2},_Acc)->
                  erase_prime(M,Key),
                  _Acc;
             (_,_Acc)->
                  _Acc
          end,[],M).
is_uncov_zero(M, CRow,CCol, Key={X, Y})->
    is_zero(M,Key) and
        not is_covered(CCol,Y) and
        not is_covered(CRow, X).

is_zero(M,Key)->
    is_true(M,Key, 0).
is_true(M, Key,Val)->
    case lookup(M,Key) of
        {Key, Val}->
            true;
        _ ->
            false
    end.
is_prime(M, Key)->
    is_true(M, Key, 2).
is_star(M, Key)->
    is_true(M, Key,1).
is_covered(M, Key)->
    is_true(M, Key, true).

prime_in_row(M,X,N)->
    PrimeFun= fun(Mat, Key={_X,Y})->
                        case is_prime(Mat, Key) of
                            true->
                                Y;
                            false ->
                                false
                        end
              end,
    find_in_row(PrimeFun, M,X,N).
star_in_row(M,X,N)->
    StarFun= fun(Mat, Key={_X,Y})->
                     case is_star(Mat, Key) of
                         true->
                             Y;
                         false ->
                             false
                     end
             end,
    find_in_row(StarFun, M,X,N).
star_in_col(M,Y,N)->
    StarFun= fun(Mat, Key={X,_Y})->
                     case is_star(Mat, Key) of
                         true->
                             X;
                         false ->
                             false
                     end
             end,
    find_in_col(StarFun, M,Y,N).

find_in_row(Fun,M,X,N)->
    find_in_row_h(Fun,M, {X, 0},N).

find_in_row_h(_Fun,_M, _Key={_X,N},N)->
    false;
find_in_row_h(Fun,M, Key={X, Y},N) ->
    case Fun(M,Key) of
        false->
            find_in_row_h(Fun,M,{X,Y+1} ,N);
        Val ->
            Val
    end.

find_in_col(Fun,M,Y,N)->
    find_in_col_h(Fun,M, {0,Y},N).

find_in_col_h(_Fun,_M, _Key={N,_Y},N)->
    false;
find_in_col_h(Fun,M, Key={X, Y},N) ->
    case Fun(M,Key) of
        false->
            find_in_col_h(Fun,M,{X+1,Y} ,N);
        Val ->
            Val
    end.
zero_in_rows(M, CRow, CCol, N)->
    zero_in_rows_h(M, CRow, CCol, 0, N).
zero_in_rows_h(_M, _CRow, _CCol, N, N) ->
    false;
zero_in_rows_h(M, CRow, CCol, X, N) ->
    case zero_in_row(M, CRow, CCol, X, N) of
        false->
            zero_in_rows_h(M, CRow, CCol, X+1, N);
        Y->
            {X,Y}
    end.
zero_in_row(M, CRow, CCol, X, N)->
    zero_in_row_h(M, CRow, CCol ,{X ,0}, N).
zero_in_row_h(_M, _CRow, _CCol ,{_X, N}, N) ->
    false;
zero_in_row_h(M, CRow, CCol, {X, Y}, N)->
    case is_uncov_zero(M, CRow,
                           CCol,{X,Y}) of
        true->
            Y;
        false->
            zero_in_row_h(M, CRow,
                          CCol, {X,Y+1},N)
    end.

covered_cols(StarMat, CCol)->
    CoverFun = fun({{_X,Y},1}, Acc)->
                       case is_covered(CCol,Y) of
                           false->
                               cover(CCol, Y),
                               Acc+1;
                           ture->
                               Acc
                       end;
                   (_,Acc)->
                       Acc
               end,
    foldl(CoverFun, 0, StarMat).


min(M)->
    MinFun =fun({_,Value},Min)when Value<Min->
                    Value;
               (_,Min)->
                    Min
            end,
    foldl(MinFun,0,M).
max(M)->
    MaxFun =fun({_,Value},Max)when Value>Max->
                    Value;
               (_,Max)->
                    Max
            end,
    foldl(MaxFun,0,M).


min_uncov(M, CRow, CCol)->
    MinUFun = fun({{X,Y}, Val},Min) when (Val /= 0)
                                         and (Val<Min)->
                      case
                          is_covered(CRow,X) and
                          is_covered(CCol,Y) of
                          false->
                              Val;
                          _->
                              Min
                      end;
                 (_,Min) ->
                      Min
              end,
    foldl(MinUFun, [], M).

sub_min(M,N)->
    sub_min_h(M,0,N).
sub_min_h(_M,N,N)->
    ok;
sub_min_h(M,X,N) ->
    ok = sub_min_row(M,get_min(M,X,N),X,N),
    sub_min_h(M,X+1,N).


sub_min_row(M, Min, X,N)->
    sub_min_row_h(M, Min, {X,0},N).
sub_min_row_h(_M, _Min, {_X,N},N)->
    ok;
sub_min_row_h(M,Min, Key={X,Y},N)->
    {Key,Val}= lookup(M, Key),
    NewVal = Val - Min,
    update(M,{Key, NewVal}),
    sub_min_row_h(M, Min, {X,Y+1},N).
get_min(M,X,N)->
    get_min_2(M,{X,0},N,[]).
get_min_2(_M,{_X,N},N,Min)->
    Min;
get_min_2(M,Key={X,Y},N,Min)->
    NewMin = case lookup(M, Key) of
                {Key, Val} when Val < Min ->
                     get_min_2(M,{X,Y+1},N,Val);
                _->
                     Min
             end,
    get_min_2(M,{X,Y+1},N,NewMin).


fold_rows(M,N,Fun,Acc)->
    fold_rows_help(M,N,{0,0},Fun,Acc).
fold_rows_help(_M,N,{N,_},_Fun,Acc)->
    Acc;
fold_rows_help(M,N,{X,N},Fun,Acc)->
    fold_rows_help(M,N,{X+1,0},Fun, Acc);
fold_rows_help(M,N,Key={X,Y},Fun,Acc) ->
    fold_rows_help(M,N,{X,Y+1},Fun, Fun(lookup(M,Key),Acc)).

star_zeros(M,N, StarMat, CRow, CCol)->
    CountFun= fun({{X,Y}, 0},Acc)->
                      case is_covered(CRow, X) or
                          is_covered(CCol, Y) of
                          false->
                              cover(CCol, Y),
                              cover(CRow, X),
                              star(StarMat,{X,Y}),
                              Acc+1;
                          _->
                              Acc
                      end;
                 (_,Acc) ->
                      Acc
              end,
    fold_rows(M,N,CountFun,0).
