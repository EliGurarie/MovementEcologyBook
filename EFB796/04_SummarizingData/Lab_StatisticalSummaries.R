## ----setup, include=FALSE---------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = TRUE, cache = TRUE, message = FALSE, warning = FALSE)
require(knitr)


## ----loadPackages, eval = TRUE, echo = TRUE, warning = FALSE, results = "hide", message=FALSE-------
# create packages list
packages <- c("plyr","dplyr","magrittr", 
              "lubridate","sf", "mapview","ggplot2")
sapply(packages, require, character = TRUE)


## ----echo = -1--------------------------------------------------------------------------------------
load("../data/elk_processed.rda")


## ----Structure--------------------------------------------------------------------------------------
str(elk_gps)


## ---------------------------------------------------------------------------------------------------
elk_df <- elk_gps |> data.frame(elk_sf |> st_transform(32611) |> st_coordinates())


## ----AddYday, eval = FALSE--------------------------------------------------------------------------
# elk_gps <- elk_gps |> mutate(doy = yday(datetime),
#                             Month = month(datetime),
#                             Year = year(datetime),
#                             Hour = hour(datetime))
# head(elk_gps)


## ----echo = FALSE-----------------------------------------------------------------------------------
elk_df <- elk_df |> mutate(doy = yday(datetime),
                            Month = month(datetime),
                            Year = year(datetime))
kable(head(elk_df))


## ----OrderTime--------------------------------------------------------------------------------------
elk_df <- elk_df |> plyr::arrange(id, datetime)


## ----IDs--------------------------------------------------------------------------------------------
length(unique(elk_df$id)) 
unique(elk_df$id)


## ---------------------------------------------------------------------------------------------------
# missing longitude
table(is.na(elk_df$lon))

# missing latitude
table(is.na(elk_df$lat))


## ----echo = 2, eval = 1-----------------------------------------------------------------------------
table(elk_df$id) |> data.frame()
table(elk_df$id)


## ----mean-------------------------------------------------------------------------------------------
mean(table(elk_df$id))


## ----sd---------------------------------------------------------------------------------------------
sd(table(elk_df$id))


## ---------------------------------------------------------------------------------------------------
min(table(elk_df$id))
max(table(elk_df$id))


## ---------------------------------------------------------------------------------------------------
summary(as.numeric(table(elk_df$id)))


## ---------------------------------------------------------------------------------------------------
hist(table(elk_df$id))


## ---------------------------------------------------------------------------------------------------
GP1 <- subset(elk_df, id =="GP1")
diff(range(GP1$datetime))


## ---------------------------------------------------------------------------------------------------
diff(range(GP1$datetime)) # auto-picks hours
difftime(max(GP1$datetime), min(GP1$datetime), units = "weeks")


## ---------------------------------------------------------------------------------------------------
diff(range(elk_df$datetime), units = "days")


## ---------------------------------------------------------------------------------------------------
(diff(range(elk_df$datetime), units = "days") |> as.numeric())/365.25


## ----eval = FALSE-----------------------------------------------------------------------------------
# elk_df |> ddply("id", plyr::summarize,
#                  start= min(datetime), end = max(datetime)) |>
#   mutate(duration = difftime(end, start, units = "days"))


## ----echo = FALSE-----------------------------------------------------------------------------------
elk_df |> ddply("id", plyr::summarize, 
                 start= min(datetime),
                 end = max(datetime)) |> 
  mutate(duration = difftime(end, start, units = "days")) |>
  head() |> kable()


## ---------------------------------------------------------------------------------------------------
elk_df |> dplyr::group_by(id) |>
  dplyr::summarize(start = min(datetime), 
                   end = max(datetime), 
                   duration = difftime(end, start, units = "days"))


## ---------------------------------------------------------------------------------------------------
elk_summary <- elk_df |> group_by(id) |> 
  dplyr::summarize(time_range = difftime(max(datetime), min(datetime), units ="days"))

elk_summary$time_range |> as.numeric() |> summary()


## ----Timerange--------------------------------------------------------------------------------------
n.summary <- elk_df |> group_by(id) |> 
  summarize(start = min(datetime), end = max(datetime)) 


## ----elkDurations-----------------------------------------------------------------------------------
require(ggplot2)
ggplot(n.summary, aes(y = id, xmin = start, xmax = end)) + 
    geom_linerange() 


## ---------------------------------------------------------------------------------------------------
n.summary <- elk_df |> group_by(id) |> 
  summarize(start = min(datetime), end = max(datetime)) %>% 
  arrange(start) |> 
  mutate(id = factor(id, levels = as.character(id)))
  
ggplot(n.summary, aes(y = id, xmin = start, xmax = end)) + 
    geom_linerange() 


## ----basePlotDurations, echo = -1-------------------------------------------------------------------
par(bty = "l", tck = 0.01, mgp = c(1,.25,0), mar = c(2,5,1,1), cex.axis = 0.8)
with(n.summary, {
  plot(start, id, xlim = range(start, end), 
       type = "n", yaxt = "n", ylab = "", xlab = "")
  segments(start, as.integer(id), end, as.integer(id), lwd = 2)
  mtext(side = 2, at = 1:nrow(n.summary), id, cex = 0.7, las = 1, line = .2)
  })


## ---------------------------------------------------------------------------------------------------
elk_days <- elk_df |> mutate(date = as.Date(datetime)) |>
  group_by(id, date) |>
  slice(1) |> arrange(id, datetime) 


## ----elkDurationPlotWithGaps------------------------------------------------------------------------
ggplot(elk_days, aes(y = id, x = datetime)) + 
    geom_point(shape = 20, size = .5, alpha = .1) 


## ----eval = FALSE-----------------------------------------------------------------------------------
# elk_df <- elk_df |> ddply("id", mutate,
#         dtime = c(NA, diff(datetime, units = "hours")))


## ----echo = FALSE-----------------------------------------------------------------------------------
head(elk_df) |> kable()


## ---------------------------------------------------------------------------------------------------
elk_df$dtime |> summary()

