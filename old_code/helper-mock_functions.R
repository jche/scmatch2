
create_mock_SL_fit <- function(){
  mock_df_to_fit = gen_one_toy()

  mock_SL_fit = get_SL_fit(df_to_fit=mock_df_to_fit,
                           X_names=c("X1","X2"),
                           Y_name="Y",
                           SL.library="SL.lm")
  return(mock_SL_fit)
}

create_mock_df_test <- function(){
  mock_df_test = gen_one_toy()
  return(mock_df_test)
}
