provide.data(stanford)
sm.survival(Age, Log.time, Status, h = 7)
sm.survival(Age, Log.time, Status, h = 7, p = 0.25,
        add = T, lty = 2)
sm.survival(Age, Log.time, Status, h = 7, p = 0.10, 
        add = T, lty = 3)
