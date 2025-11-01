# HMM Calculation Guide for Beginners

This guide explains how to manually calculate Hidden Markov Model (HMM) values for haplotype phasing, written for someone taking their first statistics course. No jargon, no Greek symbols!

## What's an HMM?

Think of the HMM like a detective solving a mystery. You're trying to figure out which of several "reference people" your DNA sample came from at each position along the chromosome.

**The Mystery:** At each DNA position, your sample matches one of the reference people perfectly, or it has an error (doesn't match any).

**The Clues:** 
- How well does the DNA match at this position?
- Which reference person did it probably match at the previous position?
- How likely is it to "switch" to a different reference person?

## Step-by-Step Example

Let's use a simple case:
- **2 reference people**: Person A and Person B
- **3 DNA positions**: Position 1, Position 2, Position 3
- **Your sample**: Has values at each position

### Initial Setup (Position 1)

**Start with equal guesses:**
- Chance it's Person A: 50% (0.5)
- Chance it's Person B: 50% (0.5)

This is called your "initial belief" - you don't know anything yet, so you guess 50-50.

### Forward Pass: Building Up Your Beliefs

**Position 1:** Sample = "T", Person A = "T", Person B = "C"

1. **Check matches:**
   - Person A matches! Score = 1.0 (perfect match)
   - Person B doesn't match. Score = 0.01 (mismatch penalty)

2. **Update beliefs:**
   - Belief for A = 0.5 × 1.0 = 0.5
   - Belief for B = 0.5 × 0.01 = 0.005

3. **Normalize (make them add to 1.0):**
   - Total = 0.5 + 0.005 = 0.505
   - Normalized A = 0.5 / 0.505 = 99%
   - Normalized B = 0.005 / 0.505 = 1%

**After Position 1:** You're 99% sure it's Person A!

**Position 2:** Sample = "?" (missing data)

1. **Transitions:** Maybe it switched to the other person?
   - Stay with same person: 95% chance
   - Switch to other person: 5% chance

2. **Calculate for each possibility:**
   - For A: (99% stayed × 95%) + (1% switched × 5%) = 94%
   - For B: (1% stayed × 95%) + (99% switched × 5%) = 6%

3. **Apply match scores:** Since data is missing, both get score = 1.0
   - Belief for A = 94% × 1.0 = 94%
   - Belief for B = 6% × 1.0 = 6%
   - (Already adds to 100%, no normalization needed)

**After Position 2:** Still 94% sure it's Person A (missing data didn't change much).

**Position 3:** Sample = "G", Person A = "G", Person B = "T"

1. **Transitions:**
   - For A: (94% stayed × 95%) + (6% switched × 5%) = 89.6%
   - For B: (6% stayed × 95%) + (94% switched × 5%) = 10.4%

2. **Apply matches:**
   - A matches: 89.6% × 1.0 = 89.6%
   - B doesn't: 10.4% × 0.01 = 0.104%

3. **Normalize:**
   - Total = 89.704%
   - A = 89.6 / 89.704 = 99.9%
   - B = 0.104 / 89.704 = 0.1%

**After Position 3:** Almost certain (99.9%) it's Person A!

### Backward Pass: Refining Earlier Guesses

The backward pass goes in reverse and uses "future information" to refine earlier beliefs.

**Start at Position 3:** (last position)
- Beta for A = 1.0
- Beta for B = 1.0
(Always start at 1.0 at the end)

**Position 2:** Use Position 3's information

For Person A:
- If it stays A → A at Pos 3: 95% × 1.0 × (A matches at 3) = 95% × 1.0 = 0.95
- If it switches to B → B at Pos 3: 5% × 1.0 × (B doesn't match at 3) = 5% × 0.01 = 0.0005
- Beta for A = 0.95 + 0.0005 = 0.9505

For Person B:
- If it stays B → B at Pos 3: 95% × 1.0 × 0.01 = 0.0095
- If it switches to A → A at Pos 3: 5% × 1.0 × 1.0 = 0.05
- Beta for B = 0.0095 + 0.05 = 0.0595

**Position 1:** Use Position 2's information

Similar calculation using the Beta values from Position 2.

### Combining Forward and Backward

Final belief = (Forward belief) × (Backward belief) / Total

This gives you the most accurate guess using both past and future information!

## Key Formulas (in Plain English)

**Forward Pass Formula:**

```
New belief = (How well this position matches) × 
             (Previous belief adjusted for possible switches)

Where:
- "How well it matches" = 1.0 if DNA matches, 0.01 if it doesn't
- "Adjusted for switches" = mostly stay with same person (95%), 
                            sometimes switch (5%)
```

**Backward Pass Formula:**

```
Beta value = Sum of all possible next states weighted by:
  - Probability of transitioning to that state
  - Beta value of that state
  - Match score at next position
```

**Transition Calculation:**

```
For each possible "next state":
  - Probability of staying = 95% (0.95)
  - Probability of switching = 5% (0.05)
  
Combined probability for state X =
  (Prob I was in X before) × (Prob stay in X) +
  (Prob I was in other states) × (Prob switch to X)
```

## How to Verify Your Calculations

**Checklist:**
1. ✅ All probabilities are between 0 and 1
2. ✅ After normalization, probabilities add to 1.0
3. ✅ Match scores are 1.0, mismatch scores are ~0.01
4. ✅ Transition probabilities (stay vs switch) add to 1.0
5. ✅ Forward pass goes left to right
6. ✅ Backward pass goes right to left
7. ✅ Final answer makes intuitive sense (matches lead to high probability)

## Common Mistakes to Avoid

1. **Forgetting to normalize** - Probabilities must add to 1.0
2. **Wrong direction** - Forward goes →, Backward goes ←
3. **Mixing up stay vs switch** - Usually ~95% stay, ~5% switch
4. **Wrong match scores** - Match = 1.0, Mismatch = 0.01 (not 0!)
5. **Missing positions** - For missing data, use score = 1.0 for everyone

## Practice Exercise

Try calculating by hand:
- 2 reference people
- 2 positions
- Position 1: Sample "A", Ref1 "A", Ref2 "T"
- Position 2: Sample "C", Ref1 "C", Ref2 "C"

What's your final belief after Position 2?

**Hint:** Both references match at Position 2, so the answer depends mainly on Position 1!

