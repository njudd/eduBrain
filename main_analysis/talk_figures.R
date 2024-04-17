#### ### ----------------- ### ####
# Nicholas Judd
# Donders Institute - LCD Lab
# 2024-04-17 https://njudd.com 


if (!require(pacman)){
  install.packages('pacman')
}

# grab functions
source("~/projects/roslaUKB/main_analysis/0_functions.R") 
# function: vec_to_fence() takes a vector and puts outliers to the fence (i.e., boxplot.stats)
# function: RDjacked() for cov corrected fuzzy RD

pacman::p_load(tidyverse, lubridate, stringr, RDHonestatchwork, ggrain)


df <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/20240219_fullset.csv")[
  !is.na(visit_date),]

# also must have some imaging data
df <-df[5 != is.na(df$CT) + is.na(df$SA) + is.na(df$WM_hyper) + is.na(df$TBV_norm) + is.na(df$CSF_norm),] # this is just a presentation it is okay there is no wFA

first_stage <- RDHonest(EduAge16 ~ running_var, data = df)

band <- first_stage$coefficients$bandwidth


df_sum <- df %>% #making sure it is imaging subjects
  group_by(running_var) %>% 
  summarise(piEdu16 = sum(EduAge16, na.rm = T)/n())

RDplt <- df_sum %>% 
  {ggplot(., aes(running_var, piEdu16)) +
      geom_point(color = "blue", alpha = .3) +
      geom_point(data=subset(., running_var > -band & running_var < band), color = "darkblue") +
      geom_point(data=subset(., running_var < -band | running_var  > band), color = "blue", alpha = .3) +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_smooth(data=subset(., running_var < 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      # geom_smooth(data=subset(., running_var < 0), method='glm',formula=y~poly(x,2),se=F, color = "darkgreen") +
      geom_smooth(data=subset(., running_var > 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = bquote('Percent staying until 16'), x = "Date of Birth in Months") + # bquote('Total Surface Area'~(mm^3))
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      ylim(c(.75, 1)) +
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())
  }



p1 <- df_sum %>% 
  {ggplot(., aes(running_var, piEdu16)) +
      geom_point() +
      labs(y = bquote('Percent staying until 16'), x = "Date of Birth in Months") + # bquote('Total Surface Area'~(mm^3))
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      ylim(c(.75, 1)) +
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())
  }


p2 <- df_sum %>% 
  {ggplot(., aes(running_var, piEdu16)) +
      geom_point(data=subset(., running_var > 0), color = "darkblue") +
      geom_point(data=subset(., running_var < 0), color = "darkgreen") +
      geom_vline(xintercept = 0,linetype="dashed") +
      labs(y = bquote('Percent staying until 16'), x = "Date of Birth in Months") + # bquote('Total Surface Area'~(mm^3))
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      ylim(c(.75, 1)) +
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())
  }


# We could just compare to the two groups... do a t-test get our soup
p2 + 
  geom_smooth(data = subset(df_sum, running_var <0), method = "lm", formula = y~1, color = "black") +
  geom_smooth(data = subset(df_sum, running_var >0), method = "lm", formula = y~1, color = "black")

# yet there is a cohort difference going on here
p2 + 
  geom_point() +
  geom_smooth(data = subset(df_sum, running_var <0), method = "lm", formula = y~x, color = "darkred") +
  geom_smooth(data = subset(df_sum, running_var <0 & running_var > -30), method = "lm", formula = y~x, color = "red")

##### could be solved my going as small as possible
# the local-randomization approach
# one month window
p1 + geom_point(data=subset(df_sum, running_var > -2 & running_var < 1), color = "red")
# five month window
p1 + geom_point(data=subset(df_sum, running_var > -6 & running_var < 5), color = "red")


# yet still has the strong assumption that those two groups are random...

# ideally we would know this function
p2 + 
  geom_point() +
  geom_smooth(data=subset(df_sum, running_var < 0), method='glm',formula=y~poly(x,4),se=F, color = "red") +
  geom_smooth(data=subset(df_sum, running_var < 0), method='glm',formula=y~poly(x,20),se=F, color = "darkred")













