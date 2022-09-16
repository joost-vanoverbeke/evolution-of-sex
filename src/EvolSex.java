

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.CombinationSampler;
import org.apache.commons.rng.sampling.distribution.MarsagliaTsangWangDiscreteSampler.Binomial;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.SharedStateDiscreteSampler;
import org.apache.commons.rng.sampling.distribution.ZigguratNormalizedGaussianSampler;
import org.apache.commons.rng.simple.RandomSource;

import java.io.*;
import java.util.Arrays;
import java.util.Comparator;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

/* class EvolvingMetacommunity
 * loops over cycles (time steps) of reproduction and dispersal
 * writes output to file */
public class EvolSex {

    static Comm comm;
    static Evol evol;
    static Run run;

    static Sites sites;

    public static void main(String[] args) throws IOException {

        System.out.println(Integer.MAX_VALUE);
        System.out.println(Long.MAX_VALUE);

        comm = new Comm();
        evol = new Evol();
        run = new Run();
        if (args.length > 0)
            Reader.readInput(args[0], comm, evol, run);

        try (PrintWriter streamOut = new PrintWriter(new FileWriter(run.fileName))) {

            long startTime = System.currentTimeMillis();
            logTitles(streamOut);

            for (int r = 0; r < run.runs; r++)
                for (int dc = 0; dc < comm.demogrCost.length; dc++)
                    for (int pc = 0; pc < comm.pChange.length; pc++)
                        for (int es = 0; es < comm.envStep.length; es++)
                            for (int dr = 0; dr < comm.dispRate.length; dr++) {

                                System.out.format("run = %d; env = %s; sex = %s; dims = %d; traits = %d; demCorr = %.2f; disp = %.4f; pChange = %.4f; step = %.4f%n",
                                        (r + 1), comm.envType, comm.sexType, comm.envDims, comm.traits, comm.demogrCost[dc], comm.dispRate[dr], comm.pChange[pc], comm.envStep[es]);

                                comm.init();
                                evol.init(comm);
                                Auxils.init(comm, evol);
                                Init init = new Init(comm);

                                sites = new Sites(comm, evol, init, dc, es, dr);

                                System.out.format("  time = %d; metacommunity N = %d; absFit = %f; relFit = %f; pSex = %f; migrCnt = %d%n",
                                        0, sites.metapopSize(), sites.absFitnessMean(), sites.relFitnessMean(), sites.pSex(), sites.migrationCounter);
                                logResults(0, streamOut, r, dc, pc, es, dr);

                                for (int t = 0; t < run.timeSteps; t++) {
                                    if (((t + 1) % (int) (1./comm.pChange[pc])) == 0)
                                        sites.changeEnvironment();
                                    sites.findMaxFitness();
                                    sites.contributionAdults();
                                    sites.reproduction();
                                    sites.disperse();

                                    if (t == 0 || ((t + 1) % run.printSteps) == 0) {
                                        System.out.format("  time = %d; metacommunity N = %d; absFit = %f; relFit = %f; pSex = %f; migrCnt = %d%n",
                                                (t + 1), sites.metapopSize(), sites.absFitnessMean(), sites.relFitnessMean(), sites.pSex(), sites.migrationCounter);
                                    }
                                    if (t == 0 || ((t + 1) % run.saveSteps) == 0) {
                                        sites.findMaxFitness();
                                        logResults(t+1, streamOut, r, dc, pc, es, dr);
                                    }
                                }
                            }

            long endTime = System.currentTimeMillis();
            System.out.println("EvolMetac took " + (endTime - startTime) +
                    " milliseconds.");
        }
    }

    static void logTitles(PrintWriter out) {
        out.print("env_type;sex_type;grid_size;patches;p_e_change;e_step;min_env;max_env;m;rho;dims;sigma_e;microsites;d;r;demogr_cost;traits;trait_loci;sex_loci;sigma_z;mu;mu_sex;omega_e;"
                + "run;time;patch;N;"
                + "p_sex_mean;p_sex_var;fitness_mean;fitness_var;abs_fitness_mean;load_mean;load_var;S_mean;S_var;"
                + "residence_distinct;residence_div;distinct_pop;div_pop");
        for (int tr = 0; tr < comm.traits; tr++)
            out.format(";dim_tr%d;e_dim_tr%d;genotype_mean_tr%d;genotype_var_tr%d;phenotype_mean_tr%d;phenotype_var_tr%d;fitness_mean_tr%d;fitness_var_tr%d;"
                            + "genotype_meta_var_tr%d;phenotype_meta_var_tr%d",
                    tr + 1, tr + 1, tr + 1, tr + 1, tr + 1, tr + 1, tr + 1, tr + 1, tr + 1, tr + 1);
        out.println("");
    }

    static void logResults(int t, PrintWriter out, int r, int dc, int pc, int es, int dr) {
        for (int p = 0; p < comm.nbrPatches; p++) {
            out.format("%s;%s;%d;%d;%f;%f;%f;%f;%f;%f;%d;%f;%d;%f;%f;%f;%d;%d;%d;%f;%f;%f;%f",
                    comm.envType, comm.sexType, comm.gridSize, comm.nbrPatches, comm.pChange[pc], comm.envStep[es], comm.minEnv, comm.maxEnv, comm.dispRate[dr], comm.rho, comm.envDims, comm.sigmaE, comm.microsites, comm.d, comm.r, comm.demogrCost[dc], comm.traits, evol.traitLoci, evol.sexLoci, evol.sigmaZ, evol.mutationRate, evol.mutationRateSex, evol.omegaE);
            out.format(";%d;%d;%d;%d",
                    r + 1, t, p + 1, sites.popSize(p));
            out.format(";"
                            + "%f;%f;%f;%f;%f;%f;%f;%f;%f;"
                            + "%f;%f;%f;%f",
                    sites.pSex(p), sites.pSexVar(p), sites.relFitnessMean(p), sites.relFitnessVar(p), sites.absFitnessMean(p), sites.relLoadMean(p), sites.relLoadVar(p), sites.selectionDiff(p), sites.selectionDiffVar(p),
                    sites.residenceDistinctMean(p),sites.residenceDivMean(p),sites.residenceDistinctPop(p),sites.residenceDivPop(p));
            for (int tr = 0; tr < comm.traits; tr++)
                out.format(";%d;%f;%f;%f;%f;%f;%f;%f;%f;%f",
                        sites.comm.traitDim[tr] + 1, sites.environment[p][sites.comm.traitDim[tr]], sites.genotypeMean(p, tr), sites.genotypeVar(p, tr), sites.phenotypeMean(p, tr), sites.phenotypeVar(p, tr), sites.traitFitnessMean(p, tr), sites.traitFitnessVar(p, tr),
                        sites.genotypeVar(tr), sites.phenotypeVar(tr));
            out.println("");
        }
    }
}


/* class Sites
 * keeps track of individuals and their attributes in microsites (microhabitats within patches)
 * implements reproduction (with inheritance and mutation) and dispersal */
class Sites {
    Comm comm;
    Evol evol;
    int totSites;

    int dcPos;
    int esPos;
    int drPos;

    int[] patch;
    boolean[] alive;
    double[][] traitPhenotype;
    double[][] traitFitness;
    double[] fitness;
    double[] pSex;

    byte[][] genotype;
    int[][] migrationGenotype;
    int migrationCounter;

    double[][] environment;
    double[] maxFitness;

    int[][] posEmpty;
    int[] nbrEmpty;
//    double[] production;
    int nbrSettled;

    int posDisp[];
    int nbrDisp;

    boolean[] sexAdults;
//    int[] endPosFathers;
//    int[][] fathersPos;
    double[] mothersProb;
    double[] fathersProb;
//    double[][] fathersCumProb;
//    double[][] mothersCumProb;


    public Sites(Comm cmm, Evol evl, Init init, int dc, int es, int dr) {
        comm = cmm;
        evol = evl;
        dcPos = dc;
        esPos = es;
        drPos = dr;

        comm.calcDispNeighbours(drPos);

        totSites = comm.nbrPatches * comm.microsites;

        patch = new int[totSites];
        alive = new boolean[totSites];
        traitPhenotype = new double[totSites][comm.traits];
        traitFitness = new double[totSites][comm.traits];
        fitness = new double[totSites];

        genotype = new byte[totSites][2 * evol.allLoci];
        migrationGenotype = new int[totSites][2 * evol.allLoci];
        migrationCounter = 2;

        pSex = new double[totSites];

        environment = new double[comm.nbrPatches][comm.envDims];
        maxFitness = new double[comm.nbrPatches];

        posEmpty = new int[comm.nbrPatches][comm.microsites];
        nbrEmpty = new int[comm.nbrPatches];
//        production = new double[comm.nbrPatches];

        posDisp = new int[totSites];

        sexAdults = new boolean[totSites];
//        endPosFathers = new int[comm.nbrPatches];
//        fathersPos = new int[comm.nbrPatches][comm.microsites];
        mothersProb = new double[totSites];
        fathersProb = new double[totSites];
//        fathersCumProb = new double[comm.nbrPatches][];
//        mothersCumProb = new double[comm.nbrPatches][totSites];

        double indGtp;
        Arrays.fill(maxFitness, 0.);

        for (int p = 0; p < comm.nbrPatches; p++) {
            if (comm.envDims >= 0) System.arraycopy(init.environment[p], 0, environment[p], 0, comm.envDims);
            for (int m = (p * comm.microsites); m < ((p + 1) * comm.microsites); m++)
                patch[m] = p;
            int[] posInds = Auxils.arraySample(init.N[p], Auxils.enumArray(p * comm.microsites, ((p + 1) * comm.microsites) - 1));
            for (int m : posInds) {
                alive[m] = true;
                fitness[m] = 1;
                for (int tr = 0; tr < comm.traits; tr++) {
                    traitFitness[m][tr] = 1;
                    indGtp = init.genotype[p][tr];
                    for (int l : evol.traitGenes[tr]) {
                        genotype[m][l] = (byte) Math.round(Auxils.random.nextDouble() * 0.5 * (Auxils.random.nextBoolean() ? -1 : 1) + indGtp);
                        migrationGenotype[m][l] = 1;
                    }

                    if (comm.sexType.equals("SWITCH")) {
                        for (int l : evol.sexGenes) {
                            genotype[m][l] = (byte) ((init.pSex < 0.5) ? 0 : 1);
                        }
                    } else {
                        for (int l : evol.sexGenes) {
                            genotype[m][l] = (byte) Math.round(Auxils.random.nextDouble() * 0.5 * (Auxils.random.nextBoolean() ? -1 : 1) + init.pSex);
                        }
                    }

                    traitPhenotype[m][tr] = calcPhenotype(m, tr);
                    traitFitness[m][tr] = calcFitness(traitPhenotype[m][tr], environment[p][comm.traitDim[tr]]);
                    fitness[m] *= traitFitness[m][tr];
                }
                if (maxFitness[p] < fitness[m])
                    maxFitness[p] = fitness[m];

                pSex[m] = Math.min(1, Math.max(0, Auxils.arrayMean(Auxils.arrayElements(genotype[m], evol.sexGenes))));
            }
        }
    }

    double calcPhenotype(int i, int tr) {
        return Auxils.arrayMean(Auxils.arrayElements(genotype[i], evol.traitGenes[tr])) + (Auxils.gaussianSampler.sample() * evol.sigmaZ);
    }

    double calcFitness(double phenot, double env) {
        return Math.exp(-(Math.pow(phenot - env, 2)) / evol.divF);
    }

//    void changeEnvironment() {
//        boolean globalEnv = comm.envType.equals("REGIONAL");
//        boolean globalChange;
//        double globalStep = 0;
//        double step;
//
//        for (int d = 0; d < comm.envDims; d++) {
//            globalChange = globalEnv && (Auxils.random.nextDouble() <= comm.pChange);
//            if (globalChange) globalStep = comm.envStep[esPos] * (Auxils.random.nextBoolean() ? -1 : 1);
//            for (int p = 0; p < comm.nbrPatches; p++) {
//                if (globalEnv ? globalChange : (Auxils.random.nextDouble() <= comm.pChange)) {
//                    step = globalEnv ? globalStep : (comm.envStep[esPos] * (Auxils.random.nextBoolean() ? -1 : 1));
//                    environment[p][d] = environment[p][d] + step;
//                    environment[p][d] = Auxils.adjustToRange(environment[p][d], comm.minEnv, comm.maxEnv);
//                    adjustFitness(p, d);
//                }
//            }
//        }
//    }

    void changeEnvironment() {
        boolean globalEnv = comm.envType.equals("REGIONAL");
        double globalStep = 0;
        double step;

        for (int d = 0; d < comm.envDims; d++) {
            globalStep = comm.envStep[esPos] * (Auxils.random.nextBoolean() ? -1 : 1);
            for (int p = 0; p < comm.nbrPatches; p++) {
                step = globalEnv ? globalStep : (comm.envStep[esPos] * (Auxils.random.nextBoolean() ? -1 : 1));
                environment[p][d] = environment[p][d] + step;
                environment[p][d] = Auxils.adjustToRange(environment[p][d], comm.minEnv, comm.maxEnv);
                adjustFitness(p, d);
            }
        }
    }

    void adjustFitness(int p, int d) {
        double oldFit;
        for (int m = (p * comm.microsites); m < ((p + 1) * comm.microsites); m++) {
            if (alive[m]) {
                oldFit = fitness[m];
                if (oldFit == 0)
                    fitness[m] = 1;
                for (int tr = 0; tr < comm.traits; tr++) {
                    if ((oldFit != 0) && (comm.traitDim[tr] == d))
                        fitness[m] /= traitFitness[m][tr];
                    if ((oldFit == 0) || (comm.traitDim[tr] == d)) {
                        traitFitness[m][tr] = calcFitness(traitPhenotype[m][tr], environment[p][comm.traitDim[tr]]);
                        fitness[m] *= traitFitness[m][tr];
                    }
                }
            }
        }
    }

    void findMaxFitness() {
        Arrays.fill(maxFitness, 0.);
        for (int i = 0; i < totSites; i++)
            if (alive[i] && maxFitness[patch[i]] < fitness[i])
                maxFitness[patch[i]] = fitness[i];
    }

    void contributionAdults() {
//        Arrays.fill(endPosFathers, 0);
        Arrays.fill(nbrEmpty, 0);

// dispersal
        nbrDisp = 0;

        double contr;
        int p;

        for (int i = 0; i < totSites; i++) {
            p = patch[i];

// dispersal
            if(Auxils.random.nextDouble() < comm.dispRate[drPos]) {
                posDisp[nbrDisp++] = i;
            }

            if (alive[i]) {
                alive[i] = Auxils.random.nextDouble() < (1 - comm.d) * (fitness[i] / maxFitness[p]);
            }
            if (alive[i]) {
                contr = 1.;
                fathersProb[i] = contr;
                sexAdults[i] = Auxils.random.nextDouble() <= pSex[i];
                if (sexAdults[i]) {
                    contr *= comm.demogrCost[dcPos];
                }
                mothersProb[i] = contr;
            } else {
                contr = 0.;
                posEmpty[p][nbrEmpty[p]++] = i;
                fathersProb[i] = contr;
                mothersProb[i] = contr;
            }
        }
    }

    void reproduction() {
        int[] posOffspring;
        double[] mothersCumProb, fathersCumProb;
        int m, f;
        double production, prod;

        for (int p = 0; p < comm.nbrPatches; p++) {
            mothersCumProb = Arrays.copyOfRange(mothersProb, (p * comm.microsites), ((p + 1) * comm.microsites));
            Auxils.arrayCumSum(mothersCumProb);
            production = mothersCumProb[comm.microsites-1] * comm.r;
            prod = (production < 1) ? ((Auxils.random.nextDouble() < production) ? 1. : 0.) : production;
            nbrSettled = Math.min(nbrEmpty[p], (int) prod);
            posOffspring = Auxils.arraySample(nbrSettled, Arrays.copyOf(posEmpty[p], nbrEmpty[p]));
            //sampling parents with replacement!
            //selfing allowed!
            for (int i = 0; i < nbrSettled; i++) {
                m = p * Auxils.randIntCumProb(mothersCumProb);
                if (sexAdults[m]) {
                    fathersCumProb = Arrays.copyOfRange(fathersProb, (p * comm.microsites), ((p + 1) * comm.microsites));
                    Auxils.arrayCumSum(fathersCumProb);
                    f = p * Auxils.randIntCumProb(fathersCumProb);
                    settle(posOffspring[i], m, f);
                } else {
                    settle(posOffspring[i], m);
                }
            }
        }
    }

    /* install newborns and inherit traits from the parent(s)
     * including mutation */
    void settle(int pos, int m, int f) {
        inherit(pos, m, f);
        mutate(pos);
        settleRest(pos);
    }

    void settle(int pos, int m) {
        inherit(pos, m);
        mutate(pos);
        settleRest(pos);
    }

     void settleRest(int pos) {
        int p = patch[pos];
        alive[pos] = true;
        fitness[pos] = 1;
        for (int tr = 0; tr < comm.traits; tr++) {
            traitPhenotype[pos][tr] = calcPhenotype(pos, tr);
            traitFitness[pos][tr] = calcFitness(traitPhenotype[pos][tr], environment[p][comm.traitDim[tr]]);
            fitness[pos] *= traitFitness[pos][tr];
        }
        if (maxFitness[p] < fitness[pos])
            maxFitness[p] = fitness[pos];
        pSex[pos] = Math.min(1, Math.max(0, Auxils.arrayMean(Auxils.arrayElements(genotype[pos], evol.sexGenes))));
    }

    /* inheritance for asexual reproduction (one parent) */
    void inherit(int posOffspring, int posParent) {
        System.arraycopy(genotype[posParent], 0, genotype[posOffspring], 0, 2 * evol.allLoci);
        System.arraycopy(migrationGenotype[posParent], 0, migrationGenotype[posOffspring], 0, 2 * evol.allLoci);
    }

    /* inheritance for sexual reproduction (two parent) */
    void inherit(int posOffspring, int posMother, int posFather) {
        for (int l = 0; l < evol.allLoci; l++) {
            if(Auxils.random.nextBoolean()) {
                genotype[posOffspring][evol.allMother[l]] = genotype[posMother][evol.allMother[l]];
                migrationGenotype[posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allMother[l]];
            } else {
                genotype[posOffspring][evol.allMother[l]] = genotype[posMother][evol.allFather[l]];
                migrationGenotype[posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allFather[l]];
            }
            if(Auxils.random.nextBoolean()) {
                genotype[posOffspring][evol.allFather[l]] = genotype[posFather][evol.allMother[l]];
                migrationGenotype[posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allMother[l]];
            } else {
                genotype[posOffspring][evol.allFather[l]] = genotype[posFather][evol.allFather[l]];
                migrationGenotype[posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allFather[l]];
            }
        }
    }

    void mutate(int posOffspring) {
        int k;
        int[] somMutLocs;

        k = Auxils.binomialSamplerSomatic.sample();
        if (k > 0) {
            CombinationSampler combinationSampler = new CombinationSampler(Auxils.random,evol.traitLoci*2, k);
            somMutLocs = Auxils.arrayElements(evol.somGenes, combinationSampler.sample());
            for (int l : somMutLocs) {
                genotype[posOffspring][l] += (Auxils.random.nextBoolean() ? -1 : 1);
            }
        }

        if (comm.sexType.equals("SWITCH")) {
            if (Auxils.random.nextDouble() <= evol.mutationRateSex) {
                for (int l : evol.sexGenes) {
                    genotype[posOffspring][l] = (byte) ((genotype[posOffspring][l] == 0) ? 1 : 0);
                }
            }
        } else {
            int[] sexMutLocs;
            double pSexTemp;
            k = Auxils.binomialSamplerSex.sample();
            if (k > 0) {
                pSexTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.sexGenes));
                CombinationSampler combinationSampler = new CombinationSampler(Auxils.random, evol.sexLoci * 2, k);
                sexMutLocs = Auxils.arrayElements(evol.sexGenes, combinationSampler.sample());
                for (int l : sexMutLocs) {
                    if (pSexTemp <= 0.)
                        genotype[posOffspring][l] += 1;
                    else if (pSexTemp >= 1.)
                        genotype[posOffspring][l] -= 1;
                    else
                        genotype[posOffspring][l] += (Auxils.random.nextBoolean() ? -1 : 1);
                }
            }
        }

    }

    // dispersal
    void disperse() {
        int oldPos, newPos;
        int[] dispShuffle = Arrays.copyOf(posDisp, nbrDisp);
        Auxils.arrayShuffle(dispShuffle);
        byte[] tempGen = Arrays.copyOf(genotype[dispShuffle[0]], 2 * evol.allLoci);
        for(int i = 1; i < nbrDisp; i++) {
            oldPos = dispShuffle[i];
            newPos = dispShuffle[i-1];
            if(alive[oldPos]) {
                System.arraycopy(genotype[oldPos], 0, genotype[newPos], 0, 2 * evol.allLoci);
                newMigrant(newPos);
                settleRest(newPos);
            }
            else {
                alive[newPos] = false;
            }
        }
        newPos = dispShuffle[nbrDisp-1];
        System.arraycopy(tempGen, 0, genotype[newPos], 0, 2 * evol.allLoci);
        newMigrant(newPos);
        settleRest(newPos);
    }

    void newMigrant(int pos) {
        Arrays.fill(migrationGenotype[pos], migrationCounter);
        migrationCounter++;
    }


    int metapopSize() {
        int tot = 0;
        for (int i = 0; i < totSites; i++)
            if (alive[i])
                tot ++;
        return tot;
    }

    int popSize(int p) {
        int tot = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                tot ++;
        return tot;
    }

    double genotypeMean(int t) {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            if (alive[i])
                mean += Auxils.arrayMean(Auxils.arrayElements(genotype[i], evol.traitGenes[t]));
        mean /= metapopSize();
        return mean;
    }

    double genotypeMean(int p, int t) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                mean += Auxils.arrayMean(Auxils.arrayElements(genotype[i], evol.traitGenes[t]));
        mean /= popSize(p);
        return mean;
    }

    double genotypeVar(int t) {
        double mean = genotypeMean(t);
        double var = 0;
        for (int i = 0; i < totSites; i++)
            if (alive[i])
                var += Math.pow(mean - Auxils.arrayMean(Auxils.arrayElements(genotype[i], evol.traitGenes[t])), 2);
        var /= metapopSize();
        return var;
    }

    double genotypeVar(int p, int t) {
        double mean = genotypeMean(p, t);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                var += Math.pow(mean - Auxils.arrayMean(Auxils.arrayElements(genotype[i], evol.traitGenes[t])), 2);
        var /= popSize(p);
        return var;
    }

    double genotypeMin(int p, int t) {
        double gtp = Auxils.arrayMean(Auxils.arrayElements(genotype[p * comm.microsites], evol.traitGenes[t]));
        double min = gtp;
        for (int i = p * comm.microsites + 1; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                gtp = Auxils.arrayMean(Auxils.arrayElements(genotype[i], evol.traitGenes[t]));
                if (gtp < min)
                    min = gtp;
            }
        }
        return min;
    }

    double genotypeMax(int p, int t) {
        double gtp = Auxils.arrayMean(Auxils.arrayElements(genotype[p * comm.microsites], evol.traitGenes[t]));
        double max = gtp;
        for (int i = p * comm.microsites + 1; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                gtp = Auxils.arrayMean(Auxils.arrayElements(genotype[i], evol.traitGenes[t]));
                if (gtp > max)
                    max = gtp;
            }
        }
        return max;
    }

    double phenotypeMean(int t) {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            if (alive[i])
                mean += traitPhenotype[i][t];
        mean /= metapopSize();
        return mean;
    }

    double phenotypeMean(int p, int t) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                mean += traitPhenotype[i][t];
        mean /= popSize(p);
        return mean;
    }

    double phenotypeVar(int t) {
        double mean = phenotypeMean(t);
        double var = 0;
        for (int i = 0; i < totSites; i++)
            if (alive[i])
                var += Math.pow(mean - traitPhenotype[i][t], 2);
        var /= metapopSize();
        return var;
    }

    double phenotypeVar(int p, int t) {
        double mean = phenotypeMean(p, t);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                var += Math.pow(mean - traitPhenotype[i][t], 2);
        var /= popSize(p);
        return var;
    }

    double phenotypeMin(int p, int t) {
        double min = traitPhenotype[p * comm.microsites][t];
        for (int i = p * comm.microsites + 1; i < (p + 1) * comm.microsites; i++) {
            if (alive[i] && traitPhenotype[i][t] < min)
                min = traitPhenotype[i][t];
        }
        return min;
    }

    double phenotypeMax(int p, int t) {
        double max = traitPhenotype[p * comm.microsites][t];
        for (int i = p * comm.microsites + 1; i < (p + 1) * comm.microsites; i++) {
            if (alive[i] && traitPhenotype[i][t] > max)
                max = traitPhenotype[i][t];
        }
        return max;
    }

    double traitFitnessMax(int p, int t) {
        double max = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i] && max < traitFitness[i][t])
                max = traitFitness[i][t];
        return max;
    }

    double traitFitnessMean(int p, int t) {
        double mean = 0;
        double max = traitFitnessMax(p, t);
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                mean += traitFitness[i][t] / max;
        mean /= popSize(p);
        return mean;
    }

    double traitFitnessVar(int p, int t) {
        double mean = traitFitnessMean(p, t);
        double max = traitFitnessMax(p, t);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                var += Math.pow(mean - traitFitness[i][t] / max, 2);
        var /= popSize(p);
        return var;
    }

    double traitFitnessMean(int p) {
        double mean = 0;
        for (int t = 0; t < comm.traits; t++)
            mean += traitFitnessMean(p, t);
        mean /= comm.traits;
        return mean;
    }

    double traitFitnessVar(int p) {
        double mean = traitFitnessMean(p);
        double var = 0;
        for (int t = 0; t < comm.traits; t++)
            var += Math.pow(mean - traitFitnessMean(p, t), 2);
        var /= popSize(p);
        return var;
    }

    double absFitnessMean() {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            if (alive[i])
                mean += fitness[i];
        mean /= metapopSize();
        return mean;
    }

    double absFitnessMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                mean += fitness[i];
        mean /= popSize(p);
        return mean;
    }

    double absFitnessVar(int p) {
        double mean = absFitnessMean(p);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                var += Math.pow(mean - fitness[i], 2);
        var /= popSize(p);
        return var;
    }

    double residenceDistinctMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                mean += Auxils.countDistinct(Auxils.arrayElements(migrationGenotype[i], evol.somGenes));
        mean /= popSize(p);
        return mean;
    }

    double residenceDistinctPop(int p) {
        int[] allGens = new int[comm.microsites*evol.traitLoci*2];
        int endGens = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                System.arraycopy(Auxils.arrayElements(migrationGenotype[i], evol.somGenes), 0, allGens, endGens, 2 * evol.traitLoci);
                endGens += 2 * evol.traitLoci;
            }
        }
        double nDistinct;
        nDistinct = Auxils.countDistinct(Arrays.copyOf(allGens, endGens));
        return nDistinct;
    }

    double residenceDivMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                mean += Auxils.divDistinct(Auxils.arrayElements(migrationGenotype[i], evol.somGenes));
        mean /= popSize(p);
        return mean;
    }

    double residenceDivPop(int p) {
        int[] allGens = new int[comm.microsites*evol.traitLoci*2];
        int endGens = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                System.arraycopy(Auxils.arrayElements(migrationGenotype[i], evol.somGenes), 0, allGens, endGens, 2 * evol.traitLoci);
                endGens += 2 * evol.traitLoci;
            }
        }
        double div;
        div = Auxils.divDistinct(Arrays.copyOf(allGens, endGens));
        return div;
    }

    double relFitnessMean() {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            if (alive[i])
                mean += (maxFitness[patch[i]] == 0) ? 0 : (fitness[i] / maxFitness[patch[i]]);
        mean /= metapopSize();
        return mean;
    }

    double relFitnessMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                mean += (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
        mean /= popSize(p);
        return mean;
    }

    double relFitnessGeom(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                mean += Math.log((maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]));
        mean /= popSize(p);
        return Math.exp(mean);
    }

    double relFitnessVar(int p) {
        double mean = relFitnessMean(p);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                var += (maxFitness[p] == 0) ? 0 : Math.pow(mean - fitness[i] / maxFitness[p], 2);
        var /= popSize(p);
        return var;
    }

    double relFitnessMin(int p) {
        double relFit, min = 1;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                relFit = (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
                if (relFit < min)
                    min = relFit;
            }
        }
        return min;
    }

    double relFitnessMax(int p) {
        double relFit, max = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                relFit = (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
                if (relFit > max)
                    max = relFit;
            }
        }
        return max;
    }

    double relLoadMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                mean += (maxFitness[p] == 0) ? 0 : (1 - fitness[i] / maxFitness[p]);
        mean /= popSize(p);
        return mean;
    }

    double relLoadGeom(int p) {
        double mean = 0;
        double min = 1;
        if (maxFitness[p] == 0)
            mean = 1;
        else {
            for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
                if (alive[i])
                    if ((fitness[i] < maxFitness[p]) && ((1 - fitness[i] / maxFitness[p]) < min))
                        min = 1 - fitness[i] / maxFitness[p];
            min = Math.exp(Math.floor(Math.log(min)));
            for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
                if (alive[i])
                    mean += Math.log(1 - fitness[i] / maxFitness[p] + min);
            mean /= popSize(p);
            mean = Math.exp(mean) - min;
        }
        return mean;
    }

    double relLoadVar(int p) {
        double mean = relLoadMean(p);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                var += (maxFitness[p] == 0) ? 0 : Math.pow(mean - (1 - fitness[i] / maxFitness[p]), 2);
        var /= popSize(p);
        return var;
    }

    double selectionDiff(int p, int t) {
        double mean = 0;
        double sum = 0;
        double fitRel;
        double fitMean = relFitnessMean(p);
        double phenotpSd = Math.sqrt(phenotypeVar(p, t));
        double SDiff;

        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
            {
                fitRel = ((maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p])) / fitMean;
                mean += traitPhenotype[i][t] * fitRel;
                sum += fitRel;
            }
        mean /= sum;
        if (phenotpSd == 0)
            SDiff = 0;
        else
            SDiff = Math.abs((mean - phenotypeMean(p, t))) / phenotpSd;
        return SDiff;
    }

    double selectionDiff(int p) {
        double mean = 0;
        for (int t = 0; t < comm.traits; t++)
            mean += selectionDiff(p, t);
        mean /= comm.traits;
        return mean;
    }

    double selectionDiffVar(int p) {
        double mean = selectionDiff(p);
        double var = 0;
        for (int t = 0; t < comm.traits; t++)
            var += Math.pow(mean - selectionDiff(p, t), 2);
        var /= popSize(p);
        return var;
    }

    double pSex() {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            if (alive[i])
                mean += pSex[i];
        mean /= metapopSize();
        return mean;
    }

    double pSex(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                mean += pSex[i];
        mean /= popSize(p);
        return mean;
    }

    double pSexVar(int p) {
        double mean = pSex(p);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                var += Math.pow(mean - pSex[i], 2);
        var /= popSize(p);
        return var;
    }
}


/* Ecological parameters/variables */
class Comm {
    String envType = "REGIONAL";
    int envDims = 1;
    int traits = 2;
    double minEnv = 0.2;
    double maxEnv = 0.8;
    double sigmaE = 0.0;
    int microsites = 600;
    double r = 5;
    double d = 0.05;
    double[] demogrCost = {0.5};

    int gridSize = 2;
    int nbrPatches = gridSize * gridSize;
    double[] pChange = {0.1};
    double[] envStep = {0.01};
    double[] dispRate = {0.01};
    double rho = 1;
    String sexType = "SWITCH";
    double pSex = 0;

    int[] traitDim;

    double[][] neighbours = new double[nbrPatches][nbrPatches];
    double[][] dispNeighbours = new double[nbrPatches][nbrPatches];

    void init() {
        nbrPatches = gridSize * gridSize;
        neighbours = new double[nbrPatches][nbrPatches];
        dispNeighbours = new double[nbrPatches][nbrPatches];
        calcDistNeighbours();

        traitDim = new int[traits];
        int dim = 0;
        for (int tr = 0; tr < traits; tr++) {
            traitDim[tr] = dim++;
            if (dim == envDims)
                dim = 0;
        }
    }

    void calcDistNeighbours() {
        for (int i = 0; i < gridSize; i++)
            for (int j = 0; j < gridSize; j++)
                for (int i2 = 0; i2 < gridSize; i2++)
                    for (int j2 = 0; j2 < gridSize; j2++) {
                        double dist = Math.sqrt(Math.pow(Math.min(Math.abs(i - i2), gridSize - Math.abs(i - i2)), 2) + Math.pow(Math.min(Math.abs(j - j2), gridSize - Math.abs(j - j2)), 2));
                        neighbours[j * gridSize + i][j2 * gridSize + i2] = dist;
                    }
    }

    void calcDispNeighbours(int dr) {
        for (int i = 0; i < nbrPatches; i++) {
            for (int j = 0; j < nbrPatches; j++)
                dispNeighbours[i][j] = (i == j) ? 0 : (rho * Math.exp(-rho * neighbours[i][j]));
            double iSum = Auxils.arraySum(dispNeighbours[i]);
            for (int j = 0; j < nbrPatches; j++)
                dispNeighbours[i][j] = (i == j) ? (1 - dispRate[dr]) : (dispRate[dr] * dispNeighbours[i][j] / iSum);
        }
    }
}


/* Evolution parameters/variables */
class Evol {
    double omegaE = 0.02;
    double divF = 1;
    int traitLoci = 20;
    int lociPerTrait = traitLoci;
    int sexLoci = 10;
    int allLoci = traitLoci + sexLoci;
    double mutationRate = 1e-4;
    double mutationRateSex = 1e-5;
    double sigmaZ = 0.01;

    int[] allMother;
    int[] allFather;
    int[] allGenes;
    int[] somMother;
    int[] somFather;
    int[] somGenes;
    int[][] traitMother;
    int[][] traitFather;
    int[][] traitGenes;
    int[] sexMother;
    int[] sexFather;
    int[] sexGenes;

    int longPos = 0;

    void init(Comm comm) {
        divF = 2 * Math.pow(Math.sqrt(comm.traits) * omegaE, 2);

        allLoci = traitLoci + sexLoci;

        lociPerTrait = traitLoci / comm.traits;

        allMother = new int[allLoci];
        allFather = new int[allLoci];
        allGenes = new int[2 * allLoci];
        somMother = new int[traitLoci];
        somFather = new int[traitLoci];
        somGenes = new int[2 * traitLoci];
        traitMother = new int[comm.traits][lociPerTrait];
        traitFather = new int[comm.traits][lociPerTrait];
        traitGenes = new int[comm.traits][2 * lociPerTrait];
        sexMother = new int[sexLoci];
        sexFather = new int[sexLoci];
        sexGenes = new int[2 * sexLoci];

        /* somatic genes */
        for (int tr = 0; tr < comm.traits; tr++) {
            for (int l = 0; l < lociPerTrait; l++) {
                longPos = l + (tr * lociPerTrait);
                traitMother[tr][l] = longPos;
                traitFather[tr][l] = traitMother[tr][l] + allLoci;
                somMother[longPos] = traitMother[tr][l];
                somFather[longPos] = somMother[longPos] + allLoci;
            }
            traitGenes[tr] = Auxils.arrayConcat(traitMother[tr], traitFather[tr]);
        }
        somGenes = Auxils.arrayConcat(somMother, somFather);

        /* sex genes */
        for (int l = 0; l < sexLoci; l++) {
            sexMother[l] = l + traitLoci;
            sexFather[l] = sexMother[l] + allLoci;
        }
        sexGenes = Auxils.arrayConcat(sexMother, sexFather);

        /* all genes */
        for (int l = 0; l < allLoci; l++) {
            allMother[l] = l;
            allFather[l] = allMother[l] + allLoci;
        }
        allGenes = Auxils.arrayConcat(allMother, allFather);
    }
}


/* run parameters */
class Run {
    int runs = 1;
    int timeSteps = 10000;
    int printSteps = 100;
    int saveSteps = 1000;
    String fileName = "output_evolution_of_sex.csv";
}


/* initialize simulation run */
class Init {
    double pSex;

    double[][] environment;
    int[] N;
    double[][] genotype;

    public Init(Comm comm) {
        double dEnv;
        environment = new double[comm.nbrPatches][comm.envDims];
        N = new int[comm.nbrPatches];
        genotype = new double[comm.nbrPatches][comm.traits];

        Arrays.fill(N, comm.microsites);

        if (comm.envType.equals("REGIONAL")) {
            for (int d = 0; d < comm.envDims; d++) {
                dEnv = comm.minEnv + (Auxils.random.nextDouble() * (comm.maxEnv - comm.minEnv));
                for (int p = 0; p < comm.nbrPatches; p++)
                    environment[p][d] = dEnv;
            }
        } else {
            for (int p = 0; p < comm.nbrPatches; p++) {
                for (int d = 0; d < comm.envDims; d++) {
                    environment[p][d] = comm.minEnv + (Auxils.random.nextDouble() * (comm.maxEnv - comm.minEnv));
                }
            }
        }

        for (int p = 0; p < comm.nbrPatches; p++) {
            for (int tr = 0; tr < comm.traits; tr++) {
                genotype[p][tr] = environment[p][comm.traitDim[tr]];
            }
        }
        pSex = comm.pSex;
    }
}


/* reading in parameter values from input file */
class Reader {
    static void readInput(String fileName, Comm comm, Evol evol, Run run) throws IOException {
        try (BufferedReader input = new BufferedReader(new FileReader(fileName))) {
            String line;
            String[] words;
            int size;
            while ((line = input.readLine()) != null) {
                words = line.trim().split("\\s+");
                switch (words[0]) {
                    case "ENVDIMS":
                        comm.envDims = Integer.parseInt(words[1]);
                        break;
                    case "TRAITS":
                        comm.traits = Integer.parseInt(words[1]);
                        break;
                    case "MINENV":
                        comm.minEnv = Double.parseDouble(words[1]);
                        break;
                    case "MAXENV":
                        comm.maxEnv = Double.parseDouble(words[1]);
                        break;
                    case "SIGMAE":
                        comm.sigmaE = Double.parseDouble(words[1]);
                        break;
                    case "MICROSITES":
                        comm.microsites = Integer.parseInt(words[1]);
                        break;
                    case "D":
                        comm.d = Double.parseDouble(words[1]);
                        break;
                    case "R":
                        comm.r = Double.parseDouble(words[1]);
                        break;
                    case "PSEX":
                        comm.pSex = Double.parseDouble(words[1]);
                        break;
                    case "COST":
                        size = Integer.parseInt(words[1]);
                        comm.demogrCost = new double[size];
                        for (int i = 0; i < size; i++)
                            comm.demogrCost[i] = Double.parseDouble(words[2 + i]);
                        break;
                    case "GRIDSIZE":
                        comm.gridSize = Integer.parseInt(words[1]);
                        break;
                    case "ENVTYPE":
                        comm.envType = words[1];
                        break;
//                    case "PCHANGE":
//                        comm.pChange = Double.parseDouble(words[1]);
//                        break;
                    case "PCHANGE":
                        size = Integer.parseInt(words[1]);
                        comm.pChange = new double[size];
                        for (int i = 0; i < size; i++)
                            comm.pChange[i] = Double.parseDouble(words[2 + i]);
                        break;
                    case "ENVSTEP":
                        size = Integer.parseInt(words[1]);
                        comm.envStep = new double[size];
                        for (int i = 0; i < size; i++)
                            comm.envStep[i] = Double.parseDouble(words[2 + i]);
                        break;
                    case "M":
                        size = Integer.parseInt(words[1]);
                        comm.dispRate = new double[size];
                        for (int i = 0; i < size; i++)
                            comm.dispRate[i] = Double.parseDouble(words[2 + i]);
                        break;
                    case "RHO":
                        comm.rho = Double.parseDouble(words[1]);
                        break;

                    case "OMEGAE":
                        evol.omegaE = Double.parseDouble(words[1]);
                        break;
                    case "TRAITLOCI":
                        evol.traitLoci = Integer.parseInt(words[1]);
                        break;
                    case "MU":
                        evol.mutationRate = Double.parseDouble(words[1]);
                        break;
                    case "SEXTYPE":
                        comm.sexType = words[1];
                        break;
                    case "SEXLOCI":
                        evol.sexLoci = Integer.parseInt(words[1]);
                        break;
                    case "MUSEX":
                        evol.mutationRateSex = Double.parseDouble(words[1]);
                        break;
                    case "SIGMAZ":
                        evol.sigmaZ = Double.parseDouble(words[1]);
                        break;

                    case "RUNS":
                        run.runs = Integer.parseInt(words[1]);
                        break;
                    case "TIMESTEPS":
                        run.timeSteps = Integer.parseInt(words[1]);
                        break;
                    case "PRINTSTEPS":
                        run.printSteps = Integer.parseInt(words[1]);
                        break;
                    case "SAVESTEPS":
                        run.saveSteps = Integer.parseInt(words[1]);
                        break;
                    case "OUTPUT":
                        run.fileName = words[1];
                        break;
                }
            }
        }
    }
}


/* Auxiliary functions for array calculations */
class Auxils {
//    static UniformRandomProvider random = RandomSource.create(RandomSource.MT_64);
//    static UniformRandomProvider random = RandomSource.create(RandomSource.JSF_64);
//    static UniformRandomProvider random = RandomSource.create(RandomSource.MSWS);
    static UniformRandomProvider random = RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP);

    static NormalizedGaussianSampler gaussianSampler = ZigguratNormalizedGaussianSampler.of(random);
    static SharedStateDiscreteSampler binomialSamplerSomatic;
    static SharedStateDiscreteSampler binomialSamplerSex;

    static void init(Comm comm, Evol evol) {
        binomialSamplerSomatic = Binomial.of(random, evol.traitLoci*2, evol.mutationRate);
        binomialSamplerSex = Binomial.of(random, evol.sexLoci*2, evol.mutationRateSex);
    }

    static void arrayShuffle(int[] array) {
        int index, temp;
        for (int i = array.length - 1; i > 0; i--) {
            index = random.nextInt(i + 1);
            temp = array[index];
            array[index] = array[i];
            array[i] = temp;
        }
    }

    static void arrayShuffle(double[] array) {
        int index;
        double temp;
        for (int i = array.length - 1; i > 0; i--) {
            index = random.nextInt(i + 1);
            temp = array[index];
            array[index] = array[i];
            array[i] = temp;
        }
    }

    static int[] arraySample(int n, int[] array) {
        int[] tempArr = array.clone();
        arrayShuffle(tempArr);
        return Arrays.copyOf(tempArr, n);
    }

    static double[] arraySample(int n, double[] array) {
        double[] tempArr = array.clone();
        arrayShuffle(tempArr);
        return Arrays.copyOf(tempArr, n);
    }

    //sampling with or without replacement
    static int[] arraySampleProb(int n, int[] array, double[] probs, boolean repl) {
        int pos;
        int[] newElements;
        int[] newArr = new int[n];
        double[] cumProbs = probs.clone();
        arrayCumSum(cumProbs);
        arrayDiv(cumProbs, cumProbs[cumProbs.length - 1]);
        for (int i = 0; i < n; i++) {
            pos = Arrays.binarySearch(cumProbs, random.nextDouble());
            pos = (pos >= 0) ? pos : (-pos - 1);
            newArr[i] = array[pos];
            if (!repl) {
                newElements = arrayConcat(enumArray(0, pos - 1), enumArray(pos + 1, array.length - 1));
                array = arrayElements(array, newElements);
                cumProbs = arrayElements(cumProbs, newElements);
                arrayCumSum(cumProbs);
                arrayDiv(cumProbs, cumProbs[cumProbs.length - 1]);
            }
        }
        return newArr;
    }

    //sampling with or without replacement
    static double[] arraySampleProb(int n, double[] array, double[] probs, boolean repl) {
        int pos;
        int[] newElements;
        double[] newArr = new double[n];
        double[] cumProbs = probs.clone();
        arrayCumSum(cumProbs);
        arrayDiv(cumProbs, cumProbs[cumProbs.length - 1]);
        for (int i = 0; i < n; i++) {
            pos = Arrays.binarySearch(cumProbs, random.nextDouble());
            pos = (pos >= 0) ? pos : (-pos - 1);
            newArr[i] = array[pos];
            if (!repl) {
                newElements = arrayConcat(enumArray(0, pos - 1), enumArray(pos + 1, array.length - 1));
                array = arrayElements(array, newElements);
                cumProbs = arrayElements(cumProbs, newElements);
                arrayCumSum(cumProbs);
                arrayDiv(cumProbs, cumProbs[cumProbs.length - 1]);
            }
        }
        return newArr;
    }

    static int randIntProb(int end, double[] probs) {
        int val;
        double[] cumProbs = Arrays.copyOf(probs, end);
        Auxils.arrayCumSum(cumProbs);
        Auxils.arrayDiv(cumProbs, cumProbs[cumProbs.length - 1]);
        val = Arrays.binarySearch(cumProbs, random.nextDouble());
        return (val >= 0) ? val : (-val - 1);
    }

    static int randIntCumProb(double[] cumProbs) {
        int val;
        val = Arrays.binarySearch(cumProbs, random.nextDouble());
        return (val >= 0) ? val : (-val - 1);
    }

    static int countDistinct(int[] arr) {
        // First sort the array so that all
        // occurrences become consecutive
        Arrays.sort(arr);

        // Traverse the sorted array
        int n = arr.length;
        int res = 0;
        for (int i = 0; i < n; i++) {

            // Move the index ahead while
            // there are duplicates
            while (i < n - 1 &&
                    arr[i] == arr[i + 1]) {
                i++;
            }
            res++;
        }
        return res;
    }

//    static int countDistinct2(int arr[]) {
//        int res = 1;
//
//        // Pick all elements one by one
//        for (int i = 1; i < arr.length; i++)
//        {
//            int j = 0;
//            for (j = 0; j < i; j++)
//                if (arr[i] == arr[j])
//                    break;
//
//            // If not printed earlier,
//            // then print it
//            if (i == j)
//                res++;
//        }
//        return res;
//    }

    static double divDistinct(int[] arr) {
        // First sort the array so that all
        // occurrences become consecutive
        Arrays.sort(arr);

        // Traverse the sorted array
        int j, n = arr.length;
        double res = 0.;
        for (int i = 0; i < n; i++) {

            // Move the index ahead while
            // there are duplicates
            j = i;
            while (i < n - 1 &&
                    arr[i] == arr[i + 1]) {
                i++;
            }
            res += Math.pow((i+1-j)/((double) n), 2.);
        }
        return 1./res;
    }

    static int[] enumArray(int from, int to) {
        int[] newArr = new int[to - from + 1];
        for (int i = 0; i < newArr.length; i++)
            newArr[i] = from++;
        return newArr;
    }

    static int[] arrayElements(int[] array, int[] pos) {
        int[] newArr = new int[pos.length];
        for (int i = 0; i < newArr.length; i++)
            newArr[i] = array[pos[i]];
        return newArr;
    }

    static byte[] arrayElements(byte[] array, int[] pos) {
        byte[] newArr = new byte[pos.length];
        for (int i = 0; i < newArr.length; i++)
            newArr[i] = array[pos[i]];
        return newArr;
    }

    static boolean[] arrayElements(boolean[] array, int[] pos) {
        boolean[] newArr = new boolean[pos.length];
        for (int i = 0; i < newArr.length; i++)
            newArr[i] = array[pos[i]];
        return newArr;
    }

    static double[] arrayElements(double[] array, int[] pos) {
        double[] newArr = new double[pos.length];
        for (int i = 0; i < newArr.length; i++)
            newArr[i] = array[pos[i]];
        return newArr;
    }

    static int[] arrayAntiElements(int[] array, int[] pos) {
        java.util.Arrays.sort(pos);
        int j = 0;
        int k = 0;
        int[] newArr = new int[array.length - pos.length];
        for (int i = 0; i < array.length; i++)
            if(i == j)
                j++;
            else
                newArr[k++] = array[i];
        return newArr;
    }

    static double arrayMean(int[] array) {
        double mean = 0;
        for (int value : array) mean += value;
        mean /= array.length;
        return mean;
    }

    static double arrayMean(byte[] array) {
        double mean = 0;
        for (byte value : array) mean += value;
        mean /= array.length;
        return mean;
    }

    static double arrayMean(boolean[] array) {
        double mean = 0;
        for (boolean b : array)
            if (b)
                mean++;
        mean /= array.length;
        return mean;
    }

    static double arrayMean(double[] array) {
        double mean = 0;
        for (double v : array) mean += v;
        mean /= array.length;
        return mean;
    }

    static double arrayMean(int[] array, int end) {
        double mean = 0;
        for (int i = 0; i < end; i++)
            mean += array[i];
        mean /= end;
        return mean;
    }

    static double arrayMean(boolean[] array, int end) {
        double mean = 0;
        for (int i = 0; i < end; i++)
            if (array[i])
                mean++;
        mean /= end;
        return mean;
    }

    static double arrayMean(double[] array, int end) {
        double mean = 0;
        for (int i = 0; i < end; i++)
            mean += array[i];
        mean /= end;
        return mean;
    }

    static int arrayMax(int[] array) {
        int max = array[0];
        for (int i = 1; i < array.length; i++)
            if (array[i] > max)
                max = array[i];
        return max;
    }

    static double arrayMax(double[] array) {
        double max = array[0];
        for (int i = 1; i < array.length; i++)
            if (array[i] > max)
                max = array[i];
        return max;
    }

    static int arrayMin(int[] array) {
        int min = array[0];
        for (int i = 1; i < array.length; i++)
            if (array[i] < min)
                min = array[i];
        return min;
    }

    static double arrayMin(double[] array) {
        double min = array[0];
        for (int i = 1; i < array.length; i++)
            if (array[i] < min)
                min = array[i];
        return min;
    }

    static int arraySum(int[] array) {
        int sum = 0;
        for (int value : array) sum += value;
        return sum;
    }

    static int arraySum(boolean[] array) {
        int sum = 0;
        for (boolean b : array)
            if (b)
                sum++;
        return sum;
    }

    static double arraySum(double[] array) {
        double sum = 0;
        for (double v : array) sum += v;
        return sum;
    }

    static int arraySum(int[] array, int end) {
        int sum = 0;
        for (int i = 0; i < end; i++)
            sum += array[i];
        return sum;
    }

    static int arraySum(boolean[] array, int end) {
        int sum = 0;
        for (int i = 0; i < end; i++)
            if (array[i])
                sum++;
        return sum;
    }

    static double arraySum(double[] array, int end) {
        double sum = 0;
        for (int i = 0; i < end; i++)
            sum += array[i];
        return sum;
    }

    static void arrayCumSum(int[] array) {
        for (int i = 1; i < array.length; i++)
            array[i] += array[i - 1];
    }

    static void arrayCumSum(double[] array) {
        for (int i = 1; i < array.length; i++)
            array[i] += array[i - 1];
    }

    static void arrayAdd(int[] array, int a) {
        for (int i = 0; i < array.length; i++)
            array[i] += a;
    }

    static void arrayAdd(double[] array, double a) {
        for (int i = 0; i < array.length; i++)
            array[i] += a;
    }

    static void arrayMult(int[] array, int a) {
        for (int i = 0; i < array.length; i++)
            array[i] *= a;
    }

    static void arrayMult(double[] array, double a) {
        for (int i = 0; i < array.length; i++)
            array[i] *= a;
    }

    static void arrayDiv(int[] array, int a) {
        for (int i = 0; i < array.length; i++)
            array[i] /= a;
    }

    static void arrayDiv(double[] array, double a) {
        for (int i = 0; i < array.length; i++)
            array[i] /= a;
    }

    static void arrayPow(double[] array, double a) {
        for (int i = 0; i < array.length; i++)
            array[i] = Math.pow(array[i], a);
    }

    static int[] arrayConcat(int[] first, int[] second) {
        int[] result = Arrays.copyOf(first, first.length + second.length);
        System.arraycopy(second, 0, result, first.length, second.length);
        return result;
    }

    static double[] arrayConcat(double[] first, double[] second) {
        double[] result = Arrays.copyOf(first, first.length + second.length);
        System.arraycopy(second, 0, result, first.length, second.length);
        return result;
    }

    static double adjustToRange(double val, double min, double max) {
        double range = max - min;
        int quot = (int) Math.floor((val - max) / range);
        double rem = mod((val - max), range);
        int minAdd = mod(quot, 2);
        int maxAdd = 1 - minAdd;
        return minAdd * (min + rem) + maxAdd * (max - rem);
    }

    static int mod(int x, int y) {
        int result = x % y;
        return result < 0 ? result + y : result;
    }

    static double mod(double x, double y) {
        double result = x % y;
        return result < 0 ? result + y : result;
    }

    static double mod(int x, double y) {
        double result = x % y;
        return result < 0 ? result + y : result;
    }

    static double mod(double x, int y) {
        double result = x % y;
        return result < 0 ? result + y : result;
    }

    /* Calculates Spearman's rank correlation coefficient, */
    static Double spearman(double [] X, double [] Y) {
        /* Error check */
        if (X == null || Y == null || X.length != Y.length) {
            return null;
        }

        /* Create Rank arrays */
        int [] rankX = getRanks(X);
        int [] rankY = getRanks(Y);

        /* Apply Spearman's formula */
        int n = X.length;
        double numerator = 0;
        for (int i = 0; i < n; i++) {
            numerator += Math.pow((rankX[i] - rankY[i]), 2);
        }
        numerator *= 6;
        return 1 - numerator / (n * ((n * n) - 1));
    }

    /* Returns a new array with ranks. Assumes unique array values. */
    public static int[] getRanks(double [] array) {
        int n = array.length;

        /* Create Pair[] and sort by values */
        Pair [] pair = new Pair[n];
        for (int i = 0; i < n; i++) {
            pair[i] = new Pair(i, array[i]);
        }
        Arrays.sort(pair, new PairValueComparator());

        /* Create and return ranks[] */
        int [] ranks = new int[n];
        int rank = 1;
        for (Pair p : pair) {
            ranks[p.index] = rank++;
        }
        return ranks;
    }

    /* A class to store 2 variables */
    public static class Pair {
        public final int index;
        public final double value;

        public Pair(int i, double v) {
            index = i;
            value = v;
        }
    }

    /* This lets us sort Pairs based on their value field */
    public static class PairValueComparator implements Comparator<Pair> {
        @Override
        public int compare(Pair p1, Pair p2) {
            if (p1.value < p2.value) {
                return -1;
            } else if (p1.value > p2.value) {
                return 1;
            } else {
                return 0;
            }
        }
    }
}


