#!/bin/bash
set -e

echo '=== SeaCache LLM Matrix Simulation ==='
echo '--- gpt2_l0_attn_qkv_sp90 ---'
./scache gpt2_l0_attn_qkv_sp90 gpt2_l0_attn_qkv_sp90 config/config.json
echo '--- gpt2_l0_attn_out_sp90 ---'
./scache gpt2_l0_attn_out_sp90 gpt2_l0_attn_out_sp90 config/config.json
echo '--- gpt2_l0_mlp_fc_sp90 ---'
./scache gpt2_l0_mlp_fc_sp90 gpt2_l0_mlp_fc_sp90 config/config.json
echo '--- gpt2_l0_mlp_proj_sp90 ---'
./scache gpt2_l0_mlp_proj_sp90 gpt2_l0_mlp_proj_sp90 config/config.json
echo '--- gpt2_l1_attn_qkv_sp90 ---'
./scache gpt2_l1_attn_qkv_sp90 gpt2_l1_attn_qkv_sp90 config/config.json
echo '--- gpt2_l1_attn_out_sp90 ---'
./scache gpt2_l1_attn_out_sp90 gpt2_l1_attn_out_sp90 config/config.json
echo '--- gpt2_l1_mlp_fc_sp90 ---'
./scache gpt2_l1_mlp_fc_sp90 gpt2_l1_mlp_fc_sp90 config/config.json
echo '--- gpt2_l1_mlp_proj_sp90 ---'
./scache gpt2_l1_mlp_proj_sp90 gpt2_l1_mlp_proj_sp90 config/config.json
echo '--- bert_l0_q_sp90 ---'
./scache bert_l0_q_sp90 bert_l0_q_sp90 config/config.json
echo '--- bert_l0_k_sp90 ---'
./scache bert_l0_k_sp90 bert_l0_k_sp90 config/config.json
echo '--- bert_l0_v_sp90 ---'
./scache bert_l0_v_sp90 bert_l0_v_sp90 config/config.json
echo '--- bert_l0_out_sp90 ---'
./scache bert_l0_out_sp90 bert_l0_out_sp90 config/config.json
echo '--- bert_l0_ffn1_sp90 ---'
./scache bert_l0_ffn1_sp90 bert_l0_ffn1_sp90 config/config.json
echo '--- bert_l0_ffn2_sp90 ---'
./scache bert_l0_ffn2_sp90 bert_l0_ffn2_sp90 config/config.json
echo '--- bert_l1_q_sp90 ---'
./scache bert_l1_q_sp90 bert_l1_q_sp90 config/config.json
echo '--- bert_l1_k_sp90 ---'
./scache bert_l1_k_sp90 bert_l1_k_sp90 config/config.json
echo '--- bert_l1_v_sp90 ---'
./scache bert_l1_v_sp90 bert_l1_v_sp90 config/config.json
echo '--- bert_l1_out_sp90 ---'
./scache bert_l1_out_sp90 bert_l1_out_sp90 config/config.json
echo '--- bert_l1_ffn1_sp90 ---'
./scache bert_l1_ffn1_sp90 bert_l1_ffn1_sp90 config/config.json
echo '--- bert_l1_ffn2_sp90 ---'
./scache bert_l1_ffn2_sp90 bert_l1_ffn2_sp90 config/config.json
echo 'Done.'
