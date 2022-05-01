sample <- c("1-234-IgG1-1_S89","10-618-IgG1-1_S98","11-278-IgG1-1_S99","12-9004-IgG1-1_S100","13-279-IgG1-1_S101","14-281-IgG1-1_S102",
            "15-291-IgG1-1_S103","16-637-IgG1-1_S104","17-641-IgG1-1_S105","18-720-IgG1-1_S106","19-257-IgG1-1_S107","2-234-IgG1-1_S90",
            "20-261-IgG1-1_S108","21-262-IgG1-1_S109","22-264-IgG1-1_S110","23-267-IgG1-1_S111","24-268-IgG1-1_S112","25-271-IgG1-1_S113",
            "26-258-IgG1-1a_S26","26-258-IgG1-1b_S70","26-258-IgG1-1c_S114","26-258-IgG1-1d_S158","26-258-IgG1-1e_S202","27-266-IgG1-1_S115",
            "28-272-IgG1-1_S116","29-273-IgG1-1_S117","3-241-IgG1-1_S91","30-275-IgG1-1_S118","31-277-IgG1-1_S119","32-234-IgG1-1_S120",
            "33-241-IgG1-1_S121","34-253-IgG1-1_S122","35-603-IgG1-1_S123","36-618-IgG1-1_S124","37-278-IgG1-1_S125","38-2004-IgG1-1_S126",
            "39-279-IgG1-1_S127","4-241-IgG1-1_S92","40-281-IgG1-1_S128","41-291-IgG1-1_S129","42-637-IgG1-1_S130","43-641-IgG1-1_S131",
            "44-720-IgG1-1_S132","5-253-IgG1-1_S93","6-253-IgG1-1_S94","7-603-IgG1-1_S95","8-603-IgG1-1_S96","9-618-IgG1-1_S97")

race <- c("Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Ayrshire","Holstein",
          "Ayrshire","Ayrshire","Ayrshire","Holstein","Ayrshire","Ayrshire",
          "Ayrshire","Ayrshire","Ayrshire","Ayrshire","Ayrshire","Ayrshire",
          "Ayrshire","Ayrshire","Holstein","Ayrshire","Ayrshire","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein")

time <- c("234","NO","278","9004","279","281",
          "291","637","641","720","NO","NO",
          "NO","NO","NO","NO","NO","NO",
          "NO","NO","NO","NO","NO","NO",
          "NO","NO","241","NO","NO","234",
          "241","253","603","618","278","9004",
          "279","NO","281","291","637","641",
          "720","253","NO","603","NO","618")

duplicate <- c("234","618","NO","NO","NO","NO",
               "NO","NO","NO","NO","NO","234",
               "NO","NO","NO","NO","NO","NO",
               "NO","NO","NO","NO","NO","NO",
               "NO","NO","241","NO","NO","NO",
               "NO","NO","NO","NO","NO","NO",
               "NO","241","NO","NO","NO","NO",
               "NO","253","253","603","603","618")

twentysix <- c("NO","NO","NO","NO","NO","NO",
               "NO","NO","NO","NO","NO","NO",
               "NO","NO","NO","NO","NO","NO",
               "TRUE","TRUE","TRUE","TRUE","TRUE","NO",
               "NO","NO","NO","NO","NO","NO",
               "NO","NO","NO","NO","NO","NO",
               "NO","NO","NO","NO","NO","NO",
               "NO","NO","NO","NO","NO","NO")
lookup <- data.frame(race = race, TS = twentysix, duplicate = duplicate, time = time)
rownames(lookup) <- sample

























sample <- c("1-234-IgG1-1_S89","10-618-IgG1-1_S98","11-278-IgG1-1_S99","12-9004-IgG1-1_S100","13-279-IgG1-1_S101","14-281-IgG1-1_S102",
            "15-291-IgG1-1_S103","16-637-IgG1-1_S104","17-641-IgG1-1_S105","18-720-IgG1-1_S106","19-257-IgG1-1_S107","2-234-IgG1-1_S90",
            "20-261-IgG1-1_S108","21-262-IgG1-1_S109","22-264-IgG1-1_S110","23-267-IgG1-1_S111","24-268-IgG1-1_S112","25-271-IgG1-1_S113",
            "26-258-IgG1-1a_S26","26-258-IgG1-1b_S70","26-258-IgG1-1c_S114","26-258-IgG1-1d_S158","26-258-IgG1-1e_S202","27-266-IgG1-1_S115",
            "28-272-IgG1-1_S116","29-273-IgG1-1_S117","3-241-IgG1-1_S91","30-275-IgG1-1_S118","31-277-IgG1-1_S119","32-234-IgG1-1_S120",
            "33-241-IgG1-1_S121","34-253-IgG1-1_S122","35-603-IgG1-1_S123","36-618-IgG1-1_S124","37-278-IgG1-1_S125","38-2004-IgG1-1_S126",
            "39-279-IgG1-1_S127","4-241-IgG1-1_S92","40-281-IgG1-1_S128","41-291-IgG1-1_S129","42-637-IgG1-1_S130","43-641-IgG1-1_S131",
            "44-720-IgG1-1_S132","5-253-IgG1-1_S93","6-253-IgG1-1_S94","7-603-IgG1-1_S95","8-603-IgG1-1_S96","9-618-IgG1-1_S97",
            "1-234-IgG2-1_S133"  , "10-618-IgG2-1_S142" , "11-278-IgG2-1_S143", 
            "12-9004-IgG2-1_S144", "13-279-IgG2-1_S145" , "14-281-IgG2-1_S146", 
            "15-291-IgG2-1_S147" , "16-637-IgG2-1_S148" , "17-641-IgG2-1_S149", 
            "18-720-IgG2-1_S150" , "19-257-IgG2-1_S151" , "2-234-IgG2-1_S134" , 
            "20-261-IgG2-1_S152" , "21-262-IgG2-1_S153" , "22-264-IgG2-1_S154", 
            "23-267-IgG2-1_S155" , "24-268-IgG2-1_S156" , "25-271-IgG2-1_S157", 
            "27-266-IgG2-1_S159" , "28-272-IgG2-1_S160" , "29-273-IgG2-1_S161", 
            "3-241-IgG2-1_S135"  , "30-275-IgG2-1_S162" , "31-277-IgG2-1_S163", 
            "32-234-IgG2-1_S164" , "33-241-IgG2-1_S165" , "34-253-IgG2-1_S166", 
            "35-603-IgG2-1_S167" , "36-618-IgG2-1_S168" , "37-278-IgG2-1_S169", 
            "38-2004-IgG2-1_S170", "39-279-IgG2-1_S171" , "4-241-IgG2-1_S136" , 
            "40-281-IgG2-1_S172" , "41-291-IgG2-1_S173" , "42-637-IgG2-1_S174", 
            "43-641-IgG2-1_S175" , "44-720-IgG2-1_S176" , "5-253-IgG2-1_S137" , 
            "6-253-IgG2-1_S138"  , "7-603-IgG2-1_S139"  , "8-603-IgG2-1_S140" , 
            "9-618-IgG2-1_S141"  , "1-234-IgG3-1_S1"    , "10-618-IgG3-1_S10" ,
            "11-278-IgG3-1_S11"  , "12-9004-IgG3-1_S12" , "13-279-IgG3-1_S13" ,
            "14-281-IgG3-1_S14"  , "15-291-IgG3-1_S15"  , "16-637-IgG3-1_S16" ,"17-641-IgG3-1_S17","18-720-IgG3-1_S18","19-257-IgG3-1_S19","2-234-IgG3-1_S2",
            "20-261-IgG3-1_S20","21-262-IgG3-1_S21","22-264-IgG3-1_S22","23-267-IgG3-1_S23","24-268-IgG3-1_S24","25-271-IgG3-1_S25",
            "27-266-IgG3-1_S27","28-272-IgG3-1_S28","29-273-IgG3-1_S29","3-241-IgG3-1_S3","30-275-IgG3-1_S30","31-277-IgG3-1_S31",
            "32-234-IgG3-1_S32","33-241-IgG3-1_S33","34-253-IgG3-1_S34","35-603-IgG3-1_S35","36-618-IgG3-1_S36","37-278-IgG3-1_S37",
            "38-2004-IgG3-1_S38","39-279-IgG3-1_S39","4-241-IgG3-1_S4","40-281-IgG3-1_S40","41-291-IgG3-1_S41","42-637-IgG3-1_S42",
            "43-641-IgG3-1_S43","44-720-IgG3-1_S44","5-253-IgG3-1_S5","6-253-IgG3-1_S6","7-603-IgG3-1_S7","8-603-IgG3-1_S8","9-618-IgG3-1_S9")

race <- c("Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Ayrshire","Holstein",
          "Ayrshire","Ayrshire","Ayrshire","Holstein","Ayrshire","Ayrshire",
          "Ayrshire","Ayrshire","Ayrshire","Ayrshire","Ayrshire","Ayrshire",
          "Ayrshire","Ayrshire","Holstein","Ayrshire","Ayrshire","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Ayrshire","Holstein",
          "Ayrshire","Ayrshire","Ayrshire","Holstein","Ayrshire","Ayrshire",
          "Ayrshire",
          "Ayrshire","Ayrshire","Holstein","Ayrshire","Ayrshire","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Ayrshire","Holstein",
          "Ayrshire","Ayrshire","Ayrshire","Holstein","Ayrshire","Ayrshire",
          "Ayrshire",
          "Ayrshire","Ayrshire","Holstein","Ayrshire","Ayrshire","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein",
          "Holstein","Holstein","Holstein","Holstein","Holstein","Holstein")

time <- c("234","NO","278","9004","279","281",
          "291","637","641","720","NO","NO",
          "NO","NO","NO","NO","NO","NO",
          "NO","NO","NO","NO","NO","NO",
          "NO","NO","241","NO","NO","234",
          "241","253","603","618","278","9004",
          "279","NO","281","291","637","641",
          "720","253","NO","603","NO","618",
          "234","NO","278","9004","279","281",
          "291","637","641","720","NO","NO",
          "NO","NO","NO","NO","NO","NO",
          "NO",
          "NO","NO","241","NO","NO","234",
          "241","253","603","618","278","9004",
          "279","NO","281","291","637","641",
          "720","253","NO","603","NO","618",
          "234","NO","278","9004","279","281",
          "291","637","641","720","NO","NO",
          "NO","NO","NO","NO","NO","NO",
          "NO",
          "NO","NO","241","NO","NO","234",
          "241","253","603","618","278","9004",
          "279","NO","281","291","637","641",
          "720","253","NO","603","NO","618")

duplicate <- c("234","618","NO","NO","NO","NO",
               "NO","NO","NO","NO","NO","234",
               "NO","NO","NO","NO","NO","NO",
               "NO","NO","NO","NO","NO","NO",
               "NO","NO","241","NO","NO","NO",
               "NO","NO","NO","NO","NO","NO",
               "NO","241","NO","NO","NO","NO",
               "NO","253","253","603","603","618",
               "234","618","NO","NO","NO","NO",
               "NO","NO","NO","NO","NO","234",
               "NO","NO","NO","NO","NO","NO",
               "NO",
               "NO","NO","241","NO","NO","NO",
               "NO","NO","NO","NO","NO","NO",
               "NO","241","NO","NO","NO","NO",
               "NO","253","253","603","603","618",
               "234","618","NO","NO","NO","NO",
               "NO","NO","NO","NO","NO","234",
               "NO","NO","NO","NO","NO","NO",
               "NO",
               "NO","NO","241","NO","NO","NO",
               "NO","NO","NO","NO","NO","NO",
               "NO","241","NO","NO","NO","NO",
               "NO","253","253","603","603","618")

Sample26_258 <- c("NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "YES","YES","YES","YES","YES","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO",
                  "NO","NO","NO","NO","NO","NO")
lookup <- data.frame(race = race, TS = Sample26_258 , duplicate = duplicate, time = time)
rownames(lookup) <- sample
