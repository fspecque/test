args <- list(
  scANVI = list(
    scANVI = c(
      # scvi.model.SCANVI
      # 'n_hidden', 'n_latent', 'n_layers', 'dropout_rate', 'dispersion',
      # 'gene_likelihood', 'linear_classifier',
      # scvi.module.SCANVAE
      'n_input', 'n_batch', 'n_labels', 'n_continuous_cov', 'n_cats_per_cov',
      'log_variational', 'y_prior', 'labels_groups', 'use_labels_groups',
      'classifier_parameters', 'use_batch_norm', 'use_layer_norm',
      # scvi.module.VAE
      'latent_distribution', 'encode_covariates', 'deeply_inject_covariates',
      'use_size_factor_key', 'use_observed_lib_size', 'library_log_means',
      'library_log_vars', 'var_activation', 'extra_encoder_kwargs',
      'extra_decoder_kwargs'#,
      # # scvi.nn.Encoder
      # 'n_output', 'n_cat_list', 'distribution', 'var_eps', 'var_activation',
      # 'return_dist',
      # # scvi.nn.FCLayers
      # 'n_in', 'n_out', 'use_activation', 'bias', 'inject_covariates', 'activation_fn'
    ),
    setup_anndata = c(
      'size_factor_key', 'categorical_covariate_keys', 'continuous_covariate_keys'
    ),
    train = c(
      'n_samples_per_label', 'check_val_every_n_epoch',
      'validation_size', 'shuffle_set_split', 'accelerator',
      'devices', 'datasplitter_kwargs', 'plan_kwargs',
      # scvi.train.Trainer
      'benchmark', 'default_root_dir', 'enable_checkpointing',
      'checkpointing_monitor', 'num_sanity_val_steps', 'enable_model_summary',
      'early_stopping', 'early_stopping_monitor', 'early_stopping_min_delta',
      'early_stopping_patience', 'early_stopping_mode', 'additional_val_metrics',
      'enable_progress_bar', 'progress_bar_refresh_rate', 'simple_progress_bar',
      'logger', 'log_every_n_steps', 'learning_rate_monitor',
      # pytorch.trainer.trainer.Trainer
      'strategy', 'num_nodes', 'precision', 'callbacks', 'fast_dev_run',
      'min_epochs', 'max_steps', 'min_steps', "max_time", 'limit_train_batches',
      'limit_val_batches', 'limit_test_batches', 'limit_predict_batches',
      'overfit_batches', 'val_check_interval', 'accumulate_grad_batches',
      'gradient_clip_val', 'gradient_clip_algorithm', 'deterministic',
      'inference_mode', 'use_distributed_sampler', 'profiler', 'detect_anomaly',
      'barebones', 'plugins', 'sync_batchnorm', 'reload_dataloaders_every_n_epochs'
    ),
    verbose = c("NOTSET"=0L, "DEBUG"=10L, "INFO"=20L, "WARNING"=30L,
                "ERROR"=40L, "CRITICAL"=50L)
  ),
  scVI = list(
    scVI = c(
      # scvi.module.VAE
      'n_input', 'n_batch', 'n_labels', 'n_continuous_cov', 'n_cats_per_cov',
      'log_variational', 'encode_covariates', 'deeply_inject_covariates',
      'use_batch_norm', 'use_layer_norm', 'use_size_factor_key',
      'use_observed_lib_size', 'library_log_means', 'library_log_vars',
      'var_activation', 'extra_encoder_kwargs',
      'extra_decoder_kwargs'
    ),
    setup_anndata = c(
      'size_factor_key', 'categorical_covariate_keys', 'continuous_covariate_keys'
    ),
    train = c(
      'accelerator', 'devices', 'validation_size', 'shuffle_set_split',
      'load_sparse_tensor', 'early_stopping', 'datasplitter_kwargs',
      'plan_kwargs', 'data_module',
      # scvi.train.Trainer
      'benchmark', 'n_samples_per_label', 'check_val_every_n_epoch',
      'default_root_dir', 'enable_checkpointing', 'checkpointing_monitor',
      'num_sanity_val_steps', 'enable_model_summary', 'early_stopping',
      'early_stopping_monitor', 'early_stopping_min_delta',
      'early_stopping_patience', 'early_stopping_mode', 'additional_val_metrics',
      'enable_progress_bar', 'progress_bar_refresh_rate', 'simple_progress_bar',
      'logger', 'log_every_n_steps', 'learning_rate_monitor',
      # pytorch.trainer.trainer.Trainer
      'strategy', 'num_nodes', 'precision', 'callbacks', 'fast_dev_run',
      'min_epochs', 'max_steps', 'min_steps', "max_time", 'limit_train_batches',
      'limit_val_batches', 'limit_test_batches', 'limit_predict_batches',
      'overfit_batches', 'val_check_interval', 'accumulate_grad_batches',
      'gradient_clip_val', 'gradient_clip_algorithm', 'deterministic',
      'inference_mode', 'use_distributed_sampler', 'profiler', 'detect_anomaly',
      'barebones', 'plugins', 'sync_batchnorm', 'reload_dataloaders_every_n_epochs'
    ),
    verbose = c("NOTSET"=0L, "DEBUG"=10L, "INFO"=20L, "WARNING"=30L,
                "ERROR"=40L, "CRITICAL"=50L)
  )
)
