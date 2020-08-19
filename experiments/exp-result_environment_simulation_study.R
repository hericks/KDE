#load(file="sim_objects.rda")

#results_sim_study <- new.env(parent=parent)

# workflow: 1. load environment
load(file="results_sim_study.rda")

# 2. set object (never uncomment already loaded objects!)
#results_sim_study$mise_kappa <- mise_kappa
#results_sim_study$mise_lambda <- mise_lambda
#results_sim_study$mise_1 <- mise_1
#results_sim_study$mise_2 <- mise_2
#results_sim_study$obj_simple_comp <- obj_simple_comp
#TODO:
#results_sim_study$plot_object_vec_1 <- plot_object_vec_1
#results_sim_study$mise_high_ns_comp <- mise_high_ns_comp
#results_sim_study$results_performance <- results_performance
#results_sim_study$plot_object_vec_2 <- plot_object_vec_2
#results_sim_study$mise_ns_comp <- mise_ns_comp


#3. save object
save(results_sim_study, file="results_sim_study.rda")


# todo load everything into global environment

parent.env(results_sim_study) <- parent.env(environment())
current <- environment()
parent.env(current) <- results_sim_study

