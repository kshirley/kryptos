################################################################################

# [1] Read in the data, basic checks:
setwd("~/misc/kryptos")
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# read in the ciphertext:
cipher_text <- readLines("ciphertext.txt")

# read in the vigenere table (the right half of the sculpture):
vigenere_table <- readLines("vigenere_table.txt")

# read in the scrabble dictionary:
scrabble <- readLines("TWL06.txt")

# read in the 2-gram frequencies:
two_gram <- fread("count_2l.txt", data.table = FALSE)
names(two_gram) <- c("letters", "frequency")

# normalize the frequencies to avoid integer overflow:
two_gram$frequency <- two_gram$frequency / min(two_gram$frequency)

# build the transition matrix:
two_mat <- matrix(0, 26, 26)
row_index <- match(substr(two_gram$letters, 1, 1), letters)
col_index <- match(substr(two_gram$letters, 2, 2), letters)
two_mat[cbind(row_index, col_index)] <- two_gram$frequency
rownames(two_mat) <- letters
colnames(two_mat) <- letters

# divide by row sums to get probabilities:
theta <- two_mat / rowSums(two_mat)
# theta

# get initial state distribution:
p <- rowSums(two_mat) / sum(two_mat)

# store the value of the vigenere table in a matrix:
v_cols <- vector("list", 26)
s <- rep(unlist(strsplit("KRYPTOSABCDEFGHIJLMNQUVWXZ", split = "")), 2)
for (i in 1:26) v_cols[[i]] <- s[1:26 + i - 1]
#for (i in 27:30) v_cols[[i]] <- s[1:26 + i - 1 - 26]
V <- do.call(cbind, v_cols)

# Let's try to decode part 1 by:
#  1. assuming it was encrypted using a vigenere cipher
#  2. trying every word in the Scrabble dictionary as the key
#  3. scoring candidate solutions using the bigram model of english

# cut out part 1 of the cipher_text:
kryptos1 <- unlist(strsplit(paste(cipher_text[1:2], collapse = ""), ""))
#key <- unlist(strsplit("PALIMPSEST", split = ""))
#key <- unlist(strsplit("ABSCISSA", split = ""))

# function to decrypt a cipher using a vigenere table and a key:
vigenere <- function(cipher, key, vigenere_table) {
  stopifnot(all(cipher %in% c(LETTERS, "?")))
  stopifnot(all(key %in% LETTERS))
  row_indices <- match(key, V[, 1])
  qmark <- cipher == "?"
  r <- row_indices[0:(length(cipher[!qmark]) - 1) %% length(row_indices) + 1]
  D <- matrix(match(V[r, ], LETTERS) - rep(match(cipher[!qmark], LETTERS), 26), 
              nrow = length(cipher[!qmark]), 
              ncol = 26)
  # one entry in each row of D equals zero, except for question marks, which are rows of NA
  out <- rep("", length(cipher))
  out[!qmark] <- V[1, (which(t(D) == 0) - 1) %% 26 + 1]
  out[qmark] <- "?"
  out
}

# test it once:
out <- vigenere(cipher = kryptos1,
                       key = unlist(strsplit("PALIMPSEST", split = "")), 
                       vigenere_table = V)

# out <- vigenere(cipher = kryptos2,
#                 key = unlist(strsplit("ABSCISSA", split = "")), 
#                 vigenere_table = V)


# function to compute the log-likelihood of a particular snippet of text:
markov_score <- function(sequence, trans_matrix = theta, initial_probs = p) {
  n <- length(sequence)
  row_indices <- match(sequence[1:(n - 1)], LETTERS)
  col_indices <- match(sequence[2:n], LETTERS)
  as.numeric(log(p[match(sequence[1], LETTERS)])) + 
    sum(log(theta[cbind(row_indices, col_indices)]), na.rm = TRUE)
}

# test it once:
# markov_score(sequence = out)

# Now, try all the words in the Scrabble dictionary for Kryptos1:
candidate <- vector("list", length(scrabble))
loglik <- numeric(length(scrabble))
system.time({
for (i in 1:length(scrabble)) {
  if (i %% 1000 == 0) print(i)
  candidate[[i]] <- vigenere(cipher = kryptos1,
                             key = unlist(strsplit(scrabble[i], split = "")), 
                             vigenere_table = V)
  loglik[i] <- markov_score(sequence = candidate[[i]])
}
})
# nice, only takes about 30 seconds

# pretty well-behaved.
hist(loglik)

# Look for the key that led to the highest likelihood solution:
scrabble[which.max(loglik)]
# [1] "PALIMPSEST"
# Great. This is, in fact, the right key.

# look at top log-likelihoods:
sort(desc(loglik))[1:20]

# look at the top 20 keys and plaintext messages.
# If we hadn't known that 'PALIMPSEST' was the correct key, could we have 
#  mistakenly concluded that another one was correct?
top_20 <- order(loglik, decreasing = TRUE)[1:20]

out <- data.frame(key = scrabble[top_20], 
                  loglik = loglik[top_20], 
                  message = sapply(candidate[top_20], paste, collapse = ""))
# OK. The second most likely message is much less likely. No ambiguity here.



# Now, let's try this for the second passage of kryptos:
kryptos2 <- unlist(strsplit(paste(cipher_text[3:14], collapse = ""), ""))

candidate <- vector("list", length(scrabble))
loglik <- numeric(length(scrabble))
system.time({
  for (i in 1:length(scrabble)) {
    if (i %% 10000 == 0) print(i)
    candidate[[i]] <- vigenere(cipher = kryptos2,
                               key = unlist(strsplit(scrabble[i], split = "")), 
                               vigenere_table = V)
    loglik[i] <- markov_score(sequence = candidate[[i]])
  }
})
# ~2 minutes for the longer message (kryptos 2)

# look at the histogram of log-likelihood values:
hist(loglik)

# Look for the key that led to the highest likelihood solution:
scrabble[which.max(loglik)]

sort(desc(loglik))[1:20]

# look at the top 20 keys and plaintext messages.
# If we hadn't known that 'ABSCISSA' was the correct key, could we have 
#  mistakenly concluded that another one was correct?
top_20 <- order(loglik, decreasing = TRUE)[1:20]

out <- data.frame(key = scrabble[top_20], 
                  loglik = loglik[top_20], 
                  message = sapply(candidate[top_20], paste, collapse = ""))

# ITWASTOTALLYINVISIBLEHOWSTHATPOSSIBLE?THEYUSEDTHEEARTHSMAGNETICFIELDXTHE
# INFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGRUUNDTOANUNKNOWNLOCATIONXDOESLANGLEY
# KNOWABOUTTHIS?THEYSHOULDITSBURIEDOUTTHERESOMEWHEREXWHOKNOWSTHEEXACTLOCATION?
# ONLYWWTHISWASHISLASTMESSAGEXTHIRTYEIGHTDEGREESFIFTYSEVENMINUTESSIXPOINTFIVE
# SECONDSNORTHSEVENTYSEVENDEGREESEIGHTMINUTESFORTYFOURSECONDSWESTIDBYROWS

# IT WAS TOTALLY INVISIBLE HOWS THAT POSSIBLE? THEY USED THE EARTHS MAGNETIC FIELD X
# THE INFORMATION WAS GATHERED AND TRANSMITTED UNDERGRUUND TO AN UNKNOWN LOCATION X 
# DOES LANGLEY KNOW ABOUT THIS? THEY SHOULD ITS BURIED OUT THERE SOMEWHERE X
# WHO KNOWS THE EXACT LOCATION? ONLY WW THIS WAS HIS LAST MESSAGE X 
# THIRTY EIGHT DEGREES FIFTY SEVEN MINUTES SIX POINT FIVE SECONDS NORTH 
# SEVENTY SEVEN DEGREES EIGHT MINUTES FORTY FOUR SECONDS WEST ID BY ROWS




### Now, try kryptos3:
kryptos3a <- unlist(strsplit(cipher_text[15:24], ""))
kryptos3b <- unlist(strsplit(cipher_text[25], ""))[1:27]
kryptos3 <- c(kryptos3a, kryptos3b)

# with a single paramters (48), set the # columns in the transposition:
n <- length(kryptos3)
o <- (seq(48, by = 48, length = n) - 1) %% n + 1
jumbled_text <- kryptos3[o]
start_index <- 4
cadence <- 4

out <- jumbled_text[(seq(start_index, by = cadence, length = n) - 1) %% n + 1]
solution <- paste(out, collapse = "")


# set a 3-d grid for (cols, start, skip)
n_cols <- 200
n_start <- 20
n_skip <- 20
trans <- expand.grid(1:n_cols, 1:n_start, 1:n_skip)
names(trans) <- c("n_cols", "start_index", "cadence")
candidate <- vector("list", nrow(trans))
ll <- numeric(nrow(trans))

system.time({
  for (i in 1:nrow(trans)) {
    if (i %% 1000 == 0) print(i)
    o <- (seq(trans$n_cols[i], by = trans$n_cols[i], length = n) - 1) %% n + 1
    candidate[[i]] <- kryptos3[o][(seq(trans$start_index[i], 
                                       by = trans$cadence[i], 
                                       length = n) - 1) %% n + 1]
    ll[i] <- markov_score(sequence = candidate[[i]])
  }
})
# takes about 10 seconds

# gather the messages into single strings instead of character vectors:
message <- sapply(candidate, paste, collapse = "")

# gather the whole solution
trans <- trans %>% mutate(ll = ll, message = message)

# look at most likely sequences:
arrange(trans, desc(ll)) %>%
  head(10)

# histogram of log-likelihoods:
hist(trans$ll)

# look at the boundary of the high likelihood cluster:
arrange(trans, desc(ll)) %>% slice(278:282)
  
# look at the combinations that led to a good solution:
filter(trans, ll > -1000) %>%
  select(-message) %>%
  arrange(n_cols, start_index, cadence) %>%
  as.data.frame()

# how many unique messages did we get using this method?
lu(trans$message[trans$ll > -1000])

sum(trans$message[trans$ll > -1000] == solution)

filter(trans, message == solution) %>%
  select(-message)
# hmmm, with n_cols = 192, start = 1, cadence = 1 you get the solution.
# seems simple! Let's unpack this.


# look at whether we actually need all three free parameters (n_cols, start_index, cadence)
m <- trans %>% group_by(message) %>%
  summarize(n = n(), 
            ind = min(start_index), 
            cad = min(cadence)) %>%
  as.data.frame()

nrow(m)
# [1] 54612
# yeah, looks like we do....


# # Look at the simplest solution:
# o <- (seq(192, by = 192, length = n) - 1) %% n + 1
# kryptos3[o]
# kryptos3[o][(seq(1, by = 1, length = n) - 1) %% n + 1]
# 
# nr <- 16
# nc <- 24
# m <- matrix("", nr, nc)
# for (i in 1:nr) {
#   m[i, ] <- kryptos3[1:nc + (i - 1)*nc]
# }
# m[is.na(m)] <- ""
# 
# 
# oo <- order(match(c("K", "R", "Y", "P", "T", "O", "S"), LETTERS))







### Part 4:
kryptos4 <- unlist(strsplit(paste(cipher_text[25:28], collapse = ""), ""))[28:124]

# Now, try all the words in the Scrabble dictionary for Kryptos4:
candidate <- vector("list", length(scrabble))
loglik <- numeric(length(scrabble))
system.time({
  for (i in 1:length(scrabble)) {
    if (i %% 1000 == 0) print(i)
    candidate[[i]] <- vigenere(cipher = kryptos4,
                               key = unlist(strsplit(scrabble[i], split = "")), 
                               vigenere_table = V)
    loglik[i] <- markov_score(sequence = candidate[[i]])
  }
})
# nice, only takes about 50 seconds

# look at the histogram of log-likelihood values:
hist(loglik)
# nothing has higher-than-random likelihood.

sort(desc(loglik))[1:20]

# Look for the key that led to the highest likelihood solution:
scrabble[which.max(loglik)]

# look at the top 20 keys and plaintext messages.
top_20 <- order(loglik, decreasing = TRUE)[1:20]

out <- data.frame(key = scrabble[top_20], 
                  loglik = loglik[top_20], 
                  message = sapply(candidate[top_20], paste, collapse = ""))
# OK, nothing here.

# check a few multiword keys:
vigenere(cipher = kryptos4,
         key = unlist(strsplit("YESWONDERFULTHINGS", split = "")), 
         vigenere_table = V)

vigenere(cipher = kryptos4,
         key = unlist(strsplit("WONDERFULTHINGS", split = "")), 
         vigenere_table = V)

vigenere(cipher = kryptos4,
         key = unlist(strsplit("SANBORN", split = "")), 
         vigenere_table = V)

vigenere(cipher = kryptos4,
         key = unlist(strsplit("JAMESSANBORN", split = "")), 
         vigenere_table = V)

vigenere(cipher = kryptos4,
         key = unlist(strsplit("JIMSANBORN", split = "")), 
         vigenere_table = V)

# Use the Berlin Clock clue:
cbind(kryptos4[64:74], c("B", "E", "R", "L", "I", "N", "C", "L", "O", "C", "K"))

# reverse engineer the vigenere key that would produce BERLIN CLOCK
col_index <- match(c("B", "E", "R", "L", "I", "N", "C", "L", "O", "C", "K"), V[1, ])
key <- rep("", 11)
for (i in 1:11) {
  key[i] <- V[which(V[, col_index[i]] == kryptos4[63 + i]), 1]
}
key
# "E" "L" "Y" "O" "I" "E" "C" "B" "A" "Q" "K"

# is the key for Berlin Clock a running key from the cipher text?
gregexpr(paste(key, collapse = ""), paste(cipher_text, collapse = ""))
# No, no match. there is no subsequence of the cipher that matches this.




# Try every transposition:
n_cols <- 90
n_start <- 20
n_skip <- 20
trans <- expand.grid(1:n_cols, 1:n_start, 1:n_skip)
names(trans) <- c("n_cols", "start_index", "cadence")
candidate <- vector("list", nrow(trans))
ll <- numeric(nrow(trans))
n4 <- length(kryptos4)

for (i in 1:nrow(trans)) {
  if (i %% 1000 == 0) print(i)
  o <- (seq(trans$n_cols[i], by = trans$n_cols[i], length = n4) - 1) %% n4 + 1
  candidate[[i]] <- kryptos4[o][(seq(trans$start_index[i], 
                                     by = trans$cadence[i], 
                                     length = n4) - 1) %% n4 + 1]
  ll[i] <- markov_score(sequence = candidate[[i]])
}

# gather the messages into single strings instead of character vectors:
message <- sapply(candidate, paste, collapse = "")

# gather the whole solution
trans <- trans %>% mutate(ll = ll, message = message)


trans %>% arrange(desc(ll)) %>%
  head(5)
hist(trans$ll)



### is the key for Berlin Clock a running key from the transposed cipher text?
search <- gregexpr(paste(key, collapse = ""), trans$message)
ss <- as.numeric(search)
table(ss)
# no matches.

# just check on the key we computed:
k <- sample(LETTERS, 97, replace = TRUE) 
k[64:74] <- c("E", "L", "Y", "O", "I", "E", "C", "B", "A", "Q", "K")
vigenere(cipher = kryptos4,
         key = k, 
         vigenere_table = V)
# well, ok. So "Berlin Clock" does show up in the right place. So yes, if
# K4 is 'finished' with the Vigenere table in the sculpture, then this 11-letter
# key would be part of the running key, or part of the keyword.





### Try a Vigenere cipher, but just look for high-likelihood letter frequencies
### i.e. not a bigram model for plaintext, but for a set of good candidates for
### a transposition to follow.

# let's try two sources of English language letter frequency:
lf_wiki <- fread("data/letter_frequencies_wikipedia.csv", data.table = FALSE)

# read in the google word frequency list:
lf_goog <- fread("data/count_1w.txt", data.table = FALSE)
names(lf_goog) <- c("word", "n")
lf_goog$len <- nchar(lf_goog$word)
x <- strsplit(lf_goog$word, split = "|")
d <- data.frame(letter = unlist(x), 
                freq = rep(lf_goog$n, lf_goog$len))
lf_goog <- d %>% group_by(letter) %>%
  summarize(freq = sum(as.numeric(freq))) %>%
  as.data.frame()
lf_goog <- mutate(lf_goog, prob = freq / sum(freq))

# join them and compare probabilities:
lf <- select(lf_wiki, letter, p_wiki = prob) %>%
  left_join(select(lf_goog, letter, p_goog = prob))

# look at the probabilities compared to each other:
lf %>% ggplot(aes(x = p_wiki, y = p_goog)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)

# read in war and peace full text:
wp <- readLines("data/book-war-and-peace.txt")
wp_vec <- unlist(strsplit(wp, split = "|"))
wp_vec <- toupper(wp_vec)
wp_vec <- wp_vec[wp_vec %in% c(letters, LETTERS)]
wp_vec <- toupper(wp_vec)
lf$p_wp <- round(as.numeric(table(wp_vec) / length(wp_vec)), 4)

lf <- rename(lf, p_google = p_goog, p_war_peace = p_wp)

# a few more comparisons:
ggplot(lf, aes(x = p_wiki, y = p_goog)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2, col = gray(0.7)) + 
  geom_text(aes(label = letter))

ggplot(lf, aes(x = p_wiki, y = p_wp)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2, col = gray(0.7)) + 
  geom_text(aes(label = letter))


# function to compute the log-likelihood of a particular snippet of text:
frequency_score <- function(sequence, letter_probs) {
  n <- length(sequence)
  sum(log(letter_probs[match(sequence, LETTERS)]), na.rm = TRUE)
}

# Now, try all the words in the Scrabble dictionary for Kryptos1:
candidate <- vector("list", length(scrabble))
loglik <- numeric(length(scrabble))
system.time({
  for (i in 1:length(scrabble)) {
    if (i %% 1000 == 0) print(i)
    candidate[[i]] <- vigenere(cipher = kryptos4,
                               key = unlist(strsplit(scrabble[i], split = "")), 
                               vigenere_table = V)
    loglik[i] <- frequency_score(sequence = candidate[[i]], letter_probs = lf$p_wp)
  }
})
# 40 seconds

hist(loglik)

# look at the quantiles of log-likelihood:
as.numeric(quantile(loglik, seq(0, 1, 0.01)))

# The question is: what should the cutoff be for where I say "there's no way 
#   this text matches the frequency distribution of English language?"

# Let's do a quick permutation test:
n_sim <- 5000
ll_sim <- numeric(n_sim)
for (i in 1:n_sim) {
  seq <- sample(LETTERS, length(kryptos4), prob = lf$p_wp, replace = TRUE)
  ll_sim[i] <- frequency_score(seq, letter_probs = lf$p_wp)
}

range(ll_sim)
# [1] -306.8252 -258.1759

range(loglik)
# [1] -440.1184 -325.7876


# wow, it looks like there is not a single Vigenere cipher solution that is within
# reasonable range of being English letter frequencies!

# Let's not bother with trying transpositions for the top-k of them.

# This seems like an effective test, though! Let's try some other cipher methods
# where we just compute the English letter frequencies, and see if any of them
# give us a sequence that could reasonably be English text.

# Let's do a quick permutation test on the text of war and peace just to check
# our null distribution:
n_sim <- 5000
ll_sim <- numeric(n_sim)
for (i in 1:n_sim) {
  seq <- wp_vec[1:97 + sample(length(wp_vec) - 98, 1)]
  ll_sim[i] <- frequency_score(seq, letter_probs = lf$p_goog)
}

range(ll_sim)
# [1] -304.9594 -261.7819


data.frame(scrabble, 
           loglik, 
           message = sapply(candidate, paste, collapse = "")) %>%
  arrange(desc(loglik)) %>%
  head(20)


tmp <- vigenere(cipher = kryptos4,
                key = unlist(strsplit("YESWONDERFULTHINGS", split = "")), 
                vigenere_table = V)
frequency_score(sequence = tmp, letter_probs = lf$p_goog)

tmp <- vigenere(cipher = kryptos4,
                key = unlist(strsplit("KRYPTOS", split = "")), 
                vigenere_table = V)
frequency_score(sequence = tmp, letter_probs = lf$p_goog)









