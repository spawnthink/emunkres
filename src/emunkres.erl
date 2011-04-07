%%% @doc
%%% This module implements the Munkres assignment algorithm (Hungarian Algorithm)
%%%
%%% @author Ahmed Omar <spawn.think@gmail.com>
-module(emunkres).
%%% @end
-export([step/2]).
-include("matrices.hrl").

step(0,[test])->
    M=?Mat:new(),
    ?Mat:from_list(M,?test_b),
    step(0,[M]);
step(0,[M])->
    N= round(math:sqrt(?Mat:get_size(M))),
    step(1,#mats{cost=M,n=N});
step(1,St) ->
    ?Mat:sub_min(?S.cost,?S.n),
    step(2,St);
step(2,St)->
    StarMat = ?Mat:new(starmat),
    CRow = ?Mat:new(cover_row),
    CCol = ?Mat:new(cover_col),
    Path = ?Mat:new(path),
    ?Mat:star_zeros(?S.cost,
                    ?S.n,
                    StarMat, CRow, CCol),
    ?Mat:uncover_all(CRow),
    ?Mat:uncover_all(CCol),
    step(3,?S{star=StarMat,
                   crow=CRow,
                   ccol=CCol,
                   path=Path});
step(3,St)->
    C=?Mat:covered_cols(?S.star,
                        ?S.ccol),
    if (C >=?S.n) ->
            ?S.star;
       true ->
            step(4,St)
    end;
step(4,St)->
    case ?Mat:zero_in_rows(?S.cost,
                           ?S.crow,
                           ?S.ccol,
                           ?S.n)
        of
        false->
            step(6,St);
        {Row,Col}->
            ?Mat:prime(?S.star,{Row,Col}),
            case ?Mat:star_in_row(?S.star,
                                  Row,
                                  ?S.n) of
                false->
                    step(5,
                         ?S{path_rc={Row,Col}});
                StarCol->
                    ?Mat:cover(?S.crow,Row),
                    ?Mat:uncover(?S.ccol,StarCol),
                    step(4,St)
            end
    end;
step(5,St=#mats{path_rc={Row,Col}})->
    ?Mat:update(?S.path,{{0,0},Row}),
    ?Mat:update(?S.path,{{0,1},Col}),
    step_5(St,0);


step(6,St)->
    MinUncov=?Mat:min_uncov(?S.cost,
                            ?S.crow,
                            ?S.ccol),
    SubFun = fun({Key={X,Y},Val},_Acc)->
                     case {?Mat:is_covered(?S.crow,X),
                           ?Mat:is_covered(?S.ccol,Y)} of
                         {true,false}->
                             ok;
                         {false,true}->
                             ok;
                         {true,_}->
                             ?Mat:update(?S.cost,{Key,Val+MinUncov});
                         {_,false} ->
                             ?Mat:update(?S.cost,{Key,Val-MinUncov})
                     end
             end,
    ?Mat:foldl(SubFun,[],?S.cost),
    step(4,St).

step_5(St,Count)->
    {_,Col} = ?Mat:lookup(?S.path,{Count,1}),
    case ?Mat:star_in_col(?S.star,
                          Col,
                         ?S.n) of
        false->
            ?Mat:convert_path(?S.star,
                              ?S.path,
                              Count+1),
            ?Mat:uncover_all(?S.crow),
            ?Mat:uncover_all(?S.ccol),
            ?Mat:erase_primes(?S.star),
            step(3,St);
        X->
            {_,V}=?Mat:lookup(?S.path,
                              {Count,1}),
            ?Mat:update(?S.path,
                        {{Count+1,0},X}),
            ?Mat:update(?S.path,
                        {{Count+1,1},V}),
            step_5(St,?Mat:check_prime(St,Count+1))
    end.



