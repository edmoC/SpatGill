library(ggplot2)
library(deSolve)

load('outputs/mvmt_pairs_h10m30.RData')

# plotting
S_H_DF['t'] <- seq(0,100)
I_H_DF['t'] <- seq(0,100)
S_M_DF['t'] <- seq(0,100)
I_M_DF['t'] <- seq(0,100)

shplot <- ggplot() + 
  theme_minimal() +
  geom_line(data=S_H_DF, aes(x=t, y=med_mh10mm1), color='lightcoral') +
  geom_line(data=S_H_DF, aes(x=t, y=med_mh7mm2), color='turquoise1') +
  geom_ribbon(data=S_H_DF, aes(x=t, ymin=bot_mh10mm1, ymax=top_mh10mm1, fill='Human Diffusion: 10\n Mosquito Diffusion: 1'), alpha=0.3) +
  geom_ribbon(data=S_H_DF, aes(x=t, ymin=bot_mh7mm2, ymax=top_mh7mm2, fill='Human Diffusion: 7\n Mosquito Diffusion: 2'), alpha=0.3) +
  labs(x='Time (days)', y='Susceptible Humans', fill='') +
  scale_color_manual(values=c('paleturquoise1', 'lightpink1'))

ihplot <- ggplot() + 
  theme_minimal() +
  geom_line(data=I_H_DF, aes(x=t, y=med_mh10mm1), color='lightcoral') +
  geom_line(data=I_H_DF, aes(x=t, y=med_mh7mm2), color='turquoise1') +
  geom_ribbon(data=I_H_DF, aes(x=t, ymin=bot_mh10mm1, ymax=top_mh10mm1, fill='Human Diffusion: 10\n Mosquito Diffusion: 1'), alpha=0.3) +
  geom_ribbon(data=I_H_DF, aes(x=t, ymin=bot_mh7mm2, ymax=top_mh7mm2, fill='Human Diffusion: 7\n Mosquito Diffusion: 2'), alpha=0.3) +
  labs(x='Time (days)', y='Infected Humans', fill='') +
  scale_color_manual(values=c('paleturquoise1', 'lightpink1'))

smplot <- ggplot() + 
  theme_minimal() +
  geom_line(data=S_M_DF, aes(x=t, y=med_mh10mm1), color='lightcoral') +
  geom_line(data=S_M_DF, aes(x=t, y=med_mh7mm2), color='turquoise1') +
  geom_ribbon(data=S_M_DF, aes(x=t, ymin=bot_mh10mm1, ymax=top_mh10mm1, fill='Human Diffusion: 10\n Mosquito Diffusion: 1'), alpha=0.3) +
  geom_ribbon(data=S_M_DF, aes(x=t, ymin=bot_mh7mm2, ymax=top_mh7mm2, fill='Human Diffusion: 7\n Mosquito Diffusion: 2'), alpha=0.3) +
  labs(x='Time (days)', y='Susceptible Mosquitoes', fill='') +
  scale_color_manual(values=c('paleturquoise1', 'lightpink1'))

implot <- ggplot() + 
  theme_minimal() +
  geom_line(data=I_M_DF, aes(x=t, y=med_mh10mm1), color='lightcoral') +
  geom_line(data=I_M_DF, aes(x=t, y=med_mh7mm2), color='turquoise1') +
  geom_ribbon(data=I_M_DF, aes(x=t, ymin=bot_mh10mm1, ymax=top_mh10mm1, fill='Human Diffusion: 10\n Mosquito Diffusion: 1'), alpha=0.3) +
  geom_ribbon(data=I_M_DF, aes(x=t, ymin=bot_mh7mm2, ymax=top_mh7mm2, fill='Human Diffusion: 7\n Mosquito Diffusion: 2'), alpha=0.3) +
  labs(x='Time (days)', y='Infected Mosquitoes', fill='') +
  scale_color_manual(values=c('paleturquoise1', 'lightpink1'))

shplot
ihplot
smplot
implot

ggsave('figures/susceptibleHumans_m13.png', shplot)
ggsave('figures/InfectedHumans_m13.png', ihplot)
ggsave('figures/susceptibleMosquitoes_m13.png', smplot)
ggsave('figures/InfectedMosquitoes_m13.png', implot)

shplot <- ggplot() + 
  theme_minimal() +
  geom_line(data=S_H_DF, aes(x=t, y=med_mh100mm10), color='lightcoral') +
  geom_line(data=S_H_DF, aes(x=t, y=med_mh70mm20), color='turquoise1') +
  geom_ribbon(data=S_H_DF, aes(x=t, ymin=bot_mh100mm10, ymax=top_mh100mm10, fill='Human Diffusion: 100\n Mosquito Diffusion: 10'), alpha=0.3) +
  geom_ribbon(data=S_H_DF, aes(x=t, ymin=bot_mh70mm20, ymax=top_mh70mm20, fill='Human Diffusion: 70\n Mosquito Diffusion: 20'), alpha=0.3) +
  labs(x='Time (days)', y='Susceptible Humans', fill='') +
  scale_color_manual(values=c('paleturquoise1', 'lightpink1'))

ihplot <- ggplot() + 
  theme_minimal() +
  geom_line(data=I_H_DF, aes(x=t, y=med_mh100mm10), color='lightcoral') +
  geom_line(data=I_H_DF, aes(x=t, y=med_mh70mm20), color='turquoise1') +
  geom_ribbon(data=I_H_DF, aes(x=t, ymin=bot_mh100mm10, ymax=top_mh100mm10, fill='Human Diffusion: 100\n Mosquito Diffusion: 10'), alpha=0.3) +
  geom_ribbon(data=I_H_DF, aes(x=t, ymin=bot_mh70mm20, ymax=top_mh70mm20, fill='Human Diffusion: 70\n Mosquito Diffusion: 20'), alpha=0.3) +
  labs(x='Time (days)', y='Infected Humans', fill='') +
  scale_color_manual(values=c('paleturquoise1', 'lightpink1'))

smplot <- ggplot() + 
  theme_minimal() +
  geom_line(data=S_M_DF, aes(x=t, y=med_mh100mm10), color='lightcoral') +
  geom_line(data=S_M_DF, aes(x=t, y=med_mh70mm20), color='turquoise1') +
  geom_ribbon(data=S_M_DF, aes(x=t, ymin=bot_mh100mm10, ymax=top_mh100mm10, fill='Human Diffusion: 100\n Mosquito Diffusion: 10'), alpha=0.3) +
  geom_ribbon(data=S_M_DF, aes(x=t, ymin=bot_mh70mm20, ymax=top_mh70mm20, fill='Human Diffusion: 70\n Mosquito Diffusion: 20'), alpha=0.3) +
  labs(x='Time (days)', y='Susceptible Mosquitoes', fill='') +
  scale_color_manual(values=c('paleturquoise1', 'lightpink1'))

implot <- ggplot() + 
  theme_minimal() +
  geom_line(data=I_M_DF, aes(x=t, y=med_mh100mm10), color='lightcoral') +
  geom_line(data=I_M_DF, aes(x=t, y=med_mh70mm20), color='turquoise1') +
  geom_ribbon(data=I_M_DF, aes(x=t, ymin=bot_mh100mm10, ymax=top_mh100mm10, fill='Human Diffusion: 100\n Mosquito Diffusion: 10'), alpha=0.3) +
  geom_ribbon(data=I_M_DF, aes(x=t, ymin=bot_mh70mm20, ymax=top_mh70mm20, fill='Human Diffusion: 70\n Mosquito Diffusion: 20'), alpha=0.3) +
  labs(x='Time (days)', y='Infected Mosquitoes', fill='') +
  scale_color_manual(values=c('paleturquoise1', 'lightpink1'))

shplot
ihplot
smplot
implot

ggsave('figures/susceptibleHumans_m130.png', shplot)
ggsave('figures/InfectedHumans_m130.png', ihplot)
ggsave('figures/susceptibleMosquitoes_m130.png', smplot)
ggsave('figures/InfectedMosquitoes_m130.png', implot)
