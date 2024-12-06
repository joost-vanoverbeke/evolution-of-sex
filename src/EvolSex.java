

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.CombinationSampler;
import org.apache.commons.rng.sampling.distribution.MarsagliaTsangWangDiscreteSampler.Binomial;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.SharedStateDiscreteSampler;
import org.apache.commons.rng.sampling.distribution.ZigguratNormalizedGaussianSampler;
import org.apache.commons.rng.simple.RandomSource;

import java.io.*;
import java.util.Arrays;

/* class EvolvingMetacommunity
 * loops over cycles (time steps) of reproduction and dispersal
 * writes output to file */
public class EvolSex {

    static Comm comm;
    static Evol evol;
    static Run run;

    static Sites sites;

    public static void main(String[] args) throws IOException {

        // System.out.println(Integer.MAX_VALUE);
        // System.out.println(Long.MAX_VALUE);

        comm = new Comm();
        evol = new Evol();
        run = new Run();
        if (args.length > 0)
            Reader.readInput(args[0], comm, evol, run);

        try (PrintWriter streamOut = new PrintWriter(new FileWriter(run.fileName), true)) {

            long startTime = System.currentTimeMillis();
            logTitles(streamOut);

            for (int r = 0; r < run.runs; r++)
            for (int dc = 0; dc < comm.demogrCost.length; dc++)
            for (int pc = 0; pc < comm.pChange.length; pc++)
            for (int es = 0; es < comm.envStep.length; es++)
            for (int dr = 0; dr < comm.dispRate.length; dr++)
            for (int ps = 0; ps < comm.pSex.length; ps++) {
                
                System.out.format("run = %d; env = %s; sex = %s; dims = %d; traits = %d; demCorr = %.2f; disp = %.4f; pChange = %.4f; step = %.4f%n",
                (r + 1), comm.envType, comm.sexType, comm.envDims, comm.traits, comm.demogrCost[dc], comm.dispRate[dr], comm.pChange[pc], comm.envStep[es]);
                
                comm.init();
                evol.init(comm);
                Auxils.init(comm, evol);
                Init init = new Init(comm, ps);
                
                sites = new Sites(comm, evol, init, dc, es, dr);
                // double[] sexSeeds = Auxils.seqArray(0., 1., 0.05);
                // for (int i = 0; i < sexSeeds.length; i++) {
                //     sexSeeds[i] = Math.round(sexSeeds[i]*100)/100.;
                // }
                // System.out.println("  probs = " + Arrays.toString(sexSeeds));
                // sites.seedSex(1, sexSeeds);
                
                System.out.format("  time = %d; metacommunity N = %d; absFit = %.2f; relFit = %.2f; pSex = %.2f; pDisp = %.5f%n",
                0, sites.metapopSize(), sites.absFitnessMean(), sites.relFitnessMean(), sites.pSex(), sites.pDisp());
                logResults(0, streamOut, r, dc, pc, es, dr, ps);
                
                for (int t = 0; t < run.timeSteps; t++) {

                    if (((t + 0) % 250) == 0)
                        // sites.seedSex(0.1, sexSeeds);
                        sites.seedSex2(0.1);

                    if (((t + 1) % (int) (1./comm.pChange[pc])) == 0)
                        sites.changeEnvironment();
                    // sites.changeEnvironment_pc(pc);
                    // sites.changeEnvironment_norm(pc);
                    sites.findMaxFitness();
                    sites.mortality();
                    sites.disperse();
                    sites.contributionAdults();
                    sites.reproduction();
                    
                    if (t == 0 || ((t + 1) % run.printSteps) == 0) {
                        System.out.format("  time = %d; metacommunity N = %d; absFit = %.2f; relFit = %.2f; pSex = %.2f; pDisp = %.5f%n",
                        (t + 1), sites.metapopSize(), sites.absFitnessMean(), sites.relFitnessMean(), sites.pSex(), sites.pDisp());
                        //                                        System.out.format("    migrationcounter = %s%n", Arrays.toString(sites.migrationCounter));
                    }
                    if (t == 0 || ((t + 1) % run.saveSteps) == 0) {
                        sites.findMaxFitness();
                        logResults(t+1, streamOut, r, dc, pc, es, dr, ps);
                    }
                }
            }

            long endTime = System.currentTimeMillis();
            System.out.println("EvolMetac took " + (endTime - startTime) +
                    " milliseconds.");
            streamOut.close();
        }
    }

    static void logTitles(PrintWriter out) {
        out.print("env_type;sex_type;init_p_sex;grid_size;patches;p_e_change;e_step;min_env;max_env;m;dims;sigma_e;microsites;d;r;demogr_cost;traits;trait_loci;sex_loci;disp_loci;sigma_z;mu;mu_sex;mu_disp;omega_e;"
                + "run;time;patch;N;"
                + "p_sex_mean;p_sex_var;p_disp_mean;p_disp_var;fitness_mean;fitness_var;abs_fitness_mean;abs_fitness_var;abs_fitness_max;load_mean;load_var;abs_contr_mean;abs_contr_var;rel_contr_var;rel_fit_var;S_mean;S_var;"
                + "distinct_pop");
        for (int tr = 0; tr < comm.traits; tr++)
            out.format(";dim_tr%d;e_dim_tr%d;genotype_mean_tr%d;genotype_var_tr%d;phenotype_mean_tr%d;phenotype_var_tr%d;"
                            + "genotype_meta_var_tr%d;phenotype_meta_var_tr%d",
                    tr + 1, tr + 1, tr + 1, tr + 1, tr + 1, tr + 1, tr + 1, tr + 1);
        out.println("");
    }

    static void logResults(int t, PrintWriter out, int r, int dc, int pc, int es, int dr, int ps) {
        for (int p = 0; p < comm.nbrPatches; p++) {
            out.format("%s;%s;%f;%d;%d;%f;%f;%f;%f;%f;%d;%f;%d;%f;%f;%f;%d;%d;%d;%d;%f;%f;%f;%f;%f",
                    comm.envType, comm.sexType, comm.pSex[ps], comm.gridSize, comm.nbrPatches, comm.pChange[pc], comm.envStep[es], comm.minEnv, comm.maxEnv, comm.dispRate[dr], comm.envDims, comm.sigmaE, comm.microsites, comm.d, comm.r, comm.demogrCost[dc], comm.traits, evol.traitLoci, evol.sexLoci, evol.dispLoci, evol.sigmaZ, evol.mutationRate, evol.mutationRateSex, evol.mutationRateDisp, evol.omegaE);
            out.format(";%d;%d;%d;%d",
                    r + 1, t, p + 1, sites.popSize(p));
            out.format(";"
                            + "%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;"
                            + "%f",
                            sites.pSex(p), sites.pSexVar(p), sites.pDisp(p), sites.pDispVar(p), sites.relFitnessMean(p), sites.relFitnessVar(p), sites.absFitnessMean(p), sites.absFitnessVar(p), sites.absFitnessMax(p), sites.relLoadMean(p), sites.relLoadVar(p), sites.absContrMean(p), sites.absContrVar(p), sites.relContrVar(p), sites.relRelFitnessVar(p), sites.selectionDiff(p), sites.selectionDiffVar(p),
                    sites.residenceDistinctPop(p));
            for (int tr = 0; tr < comm.traits; tr++)
                out.format(";%d;%f;%f;%f;%f;%f;%f;%f",
                        sites.comm.traitDim[tr] + 1, sites.environment[p][sites.comm.traitDim[tr]], sites.genotypeMean(p, tr), sites.genotypeVar(p, tr), sites.phenotypeMean(p, tr), sites.phenotypeVar(p, tr),
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
    int psPos;

    int[] patch;
    boolean[] alive;
    double[][] traitPhenotype;
    double[][] traitFitness;
    double[] fitness;
    double[] pSex;
    double[] pDisp;

    byte[][] genotype;
    int[][] migrationGenotype;
    int[] migrationCounter;

    double[][] environment;
    double[] maxFitness;

    int[] popN;

    int[][] posEmpty;
    int[] nbrEmpty;
    double[] production;
    int nbrSettled;

    int posDisp[];
    int nbrDisp;

    double[] pDispSum;

    boolean[] sexAdults;
    int[] endPosMothers;
    int[] endPosFathers;
    int[][] mothersPos;
    int[][] fathersPos;
    double[][] mothersProb;
    // double[][] fathersProb;
    double[][] mothersCumProb;
    // double[][] fathersCumProb;


    public Sites(Comm cmm, Evol evl, Init init, int dc, int es, int dr) {
        comm = cmm;
        evol = evl;
        dcPos = dc;
        esPos = es;
        drPos = dr;

        totSites = comm.nbrPatches * comm.microsites;

        patch = new int[totSites];
        alive = new boolean[totSites];
        traitPhenotype = new double[totSites][comm.traits];
        traitFitness = new double[totSites][comm.traits];
        fitness = new double[totSites];

        genotype = new byte[totSites][2 * evol.allLoci];
        migrationGenotype = new int[totSites][2 * evol.allLoci];
        migrationCounter = new int[comm.nbrPatches];
        Arrays.fill(migrationCounter, 1);

        pSex = new double[totSites];
        pDisp = new double[totSites];

        environment = new double[comm.nbrPatches][comm.envDims];
        maxFitness = new double[comm.nbrPatches];

        popN = new int[comm.nbrPatches];

        posEmpty = new int[comm.nbrPatches][comm.microsites];
        nbrEmpty = new int[comm.nbrPatches];
        production = new double[comm.nbrPatches];

        posDisp = new int[totSites];
        pDispSum = new double[comm.nbrPatches];

        sexAdults = new boolean[totSites];

        endPosMothers = new int[comm.nbrPatches];
        mothersPos = new int[comm.nbrPatches][comm.microsites];
        mothersProb = new double[comm.nbrPatches][comm.microsites];
        mothersCumProb = new double[comm.nbrPatches][comm.microsites];

        endPosFathers = new int[comm.nbrPatches];
        fathersPos = new int[comm.nbrPatches][comm.microsites];
        // fathersProb = new double[comm.nbrPatches][comm.microsites];
        // fathersCumProb = new double[comm.nbrPatches][comm.microsites];

        double indGtp;
        Arrays.fill(maxFitness, 0.);

        for (int p = 0; p < comm.nbrPatches; p++) {
            if (comm.envDims >= 0) System.arraycopy(init.environment[p], 0, environment[p], 0, comm.envDims);
            for (int m = (p * comm.microsites); m < ((p + 1) * comm.microsites); m++)
                patch[m] = p;
            int[] posInds = Auxils.arraySample(init.N[p], Auxils.enumArray(p * comm.microsites, ((p + 1) * comm.microsites) - 1));
            for (int m : posInds) {
                alive[m] = true;
                popN[p]++;
                fitness[m] = 1;
                for (int tr = 0; tr < comm.traits; tr++) {
                    traitFitness[m][tr] = 1;
                    indGtp = init.genotype[p][tr];
                    for (int l : evol.traitGenes[tr]) {
                        genotype[m][l] = (byte) Math.round(Auxils.random.nextDouble() * 0.5 * (Auxils.random.nextBoolean() ? -1 : 1) + indGtp);
                        migrationGenotype[m][l] = 0;
                    }

                    if (comm.sexType.equals("SWITCH")) {
                        for (int l : evol.sexGenes) {
                            genotype[m][l] = (byte) ((init.pSex < 0.5) ? 0 : 1);
                        }
                    } else {
                        // for (int l : evol.sexGenes) {
                        //     genotype[m][l] = (byte) Math.round(Auxils.random.nextDouble() * 0.5 * (Auxils.random.nextBoolean() ? -1 : 1) + init.pSex);
                        // }
                        for (int l = 0; l < evol.sexLoci; l++) {
                            if (l < Math.round(init.pSex*evol.sexLoci)) {
                                genotype[m][evol.sexMother[l]] = genotype[m][evol.sexFather[l]] = (byte) 1;
                            } else {
                                genotype[m][evol.sexMother[l]] = genotype[m][evol.sexFather[l]] = (byte) 0;
                            }
                        }
                    }

                    for (int l = 0; l < evol.dispLoci; l++) {
                        // if (l < Math.round(Math.pow(comm.dispRate[drPos]*evol.dispCorr, 1/evol.dispPow)*evol.dispLoci)) {
                        if (l < Math.round((1 - Math.log(comm.dispRate[drPos])/Math.log(evol.refDisp)*(1-0.2))*evol.dispLoci)) {
                                genotype[m][evol.dispMother[l]] = genotype[m][evol.dispFather[l]] = (byte) 1;
                        } else {
                            genotype[m][evol.dispMother[l]] = genotype[m][evol.dispFather[l]] = (byte) 0;
                        }
                    }

                    traitPhenotype[m][tr] = calcPhenotype(m, tr);
                    traitFitness[m][tr] = calcFitness(traitPhenotype[m][tr], environment[p][comm.traitDim[tr]]);
                    fitness[m] *= traitFitness[m][tr];
                }
                if (maxFitness[p] < fitness[m])
                    maxFitness[p] = fitness[m];
                    // pSex[m] = Math.min(1, Math.max(0, Auxils.arrayMean(Auxils.arrayElements(genotype[m], evol.sexGenes))));
                    pSex[m] = Math.min(1, Math.max(0, Auxils.arrayMean(genotype[m], evol.sexGenes)));
                    // pDisp[m] = Math.pow(Math.min(1, Math.max(0, Auxils.arrayMean(Auxils.arrayElements(genotype[m], evol.dispGenes)))),evol.dispPow)/evol.dispCorr;
                    // pDisp[m] = Math.exp(Math.log(evol.refDisp)*(1 - Math.min(1, Math.max(0, Auxils.arrayMean(Auxils.arrayElements(genotype[m], evol.dispGenes)))))/(1 - 0.2));
                    pDisp[m] = Math.exp(Math.log(evol.refDisp)*(1 - Math.min(1, Math.max(0, Auxils.arrayMean(genotype[m], evol.dispGenes))))/(1 - 0.2));

                pDispSum[p] += pDisp[m];
            }
        }
    }

    double calcPhenotype(int i, int tr) {
        return Auxils.arrayMean(genotype[i], evol.traitGenes[tr]) + (Auxils.gaussianSampler.sample() * evol.sigmaZ);
    }

    double calcFitness(double phenot, double env) {
        return Math.exp(-(Math.pow(phenot - env, 2)) / evol.divF);
    }

    void seedSex(double r, double[] probs) {
        double pS;
        for (int i = 0; i < totSites; i++) {
            if (alive[i] && Auxils.random.nextDouble() <= r) {
                // pS = Auxils.arraySample(1, Auxils.enumArray(0, 20))[0]/20.;
                // pS = Auxils.arraySample(1, probs)[0];
                pS = probs[Auxils.random.nextInt(probs.length)];
                for (int l = 0; l < evol.sexLoci; l++) {
                    if (l < Math.round(pS*evol.sexLoci)) {
                        genotype[i][evol.sexMother[l]] = genotype[i][evol.sexFather[l]] = (byte) 1;
                    } else {
                        genotype[i][evol.sexMother[l]] = genotype[i][evol.sexFather[l]] = (byte) 0;
                    }
                }
                pSex[i] = Math.min(1, Math.max(0, Auxils.arrayMean(genotype[i], evol.sexGenes)));
            }
        }
    }


    void seedSex2(double r) {
        double pS;
        for (int i = 0; i < totSites; i++) {
            if (alive[i] && Auxils.random.nextDouble() <= r) {
                // pS = Auxils.arraySample(1, Auxils.enumArray(0, 20))[0]/20.;
                // pS = Auxils.arraySample(1, probs)[0];
                pS = pSex[i] + (Auxils.random.nextBoolean() ? -1 : 1)*1./evol.sexLoci;
                pS = Math.max(Math.min(pS, 1.), 0.);
                for (int l = 0; l < evol.sexLoci; l++) {
                    if (l < Math.round(pS*evol.sexLoci)) {
                        genotype[i][evol.sexMother[l]] = genotype[i][evol.sexFather[l]] = (byte) 1;
                    } else {
                        genotype[i][evol.sexMother[l]] = genotype[i][evol.sexFather[l]] = (byte) 0;
                    }
                }
                pSex[i] = Math.min(1, Math.max(0, Auxils.arrayMean(genotype[i], evol.sexGenes)));
            }
        }
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

void changeEnvironment_pc(int pc) {
    boolean globalEnv = comm.envType.equals("REGIONAL");
    double step = 0;

    if (globalEnv) {
        for (int d = 0; d < comm.envDims; d++) {
            if (Auxils.random.nextDouble() <= comm.pChange[pc]) {
                step = comm.envStep[esPos] * (Auxils.random.nextBoolean() ? -1 : 1);
                for (int p = 0; p < comm.nbrPatches; p++) {
                    environment[p][d] += step;
                    environment[p][d] = Auxils.adjustToRange(environment[p][d], comm.minEnv, comm.maxEnv);
                    adjustFitness(p, d);
                }
            }
        }
    } else {
        for (int d = 0; d < comm.envDims; d++) {
            for (int p = 0; p < comm.nbrPatches; p++) {
                if (Auxils.random.nextDouble() <= comm.pChange[pc]) {
                    step = comm.envStep[esPos] * (Auxils.random.nextBoolean() ? -1 : 1);
                    environment[p][d] += step;
                    environment[p][d] = Auxils.adjustToRange(environment[p][d], comm.minEnv, comm.maxEnv);
                    adjustFitness(p, d);
                }
            }
        }
    }
}

void changeEnvironment_norm(int pc) {
    boolean globalEnv = comm.envType.equals("REGIONAL");
    double step = 0;

    if (globalEnv) {
        for (int d = 0; d < comm.envDims; d++) {
            if (Auxils.random.nextDouble() <= comm.pChange[pc]) {
                step = Auxils.gaussianSampler.sample() * comm.envStep[esPos];
                for (int p = 0; p < comm.nbrPatches; p++) {
                    environment[p][d] += step;
                    environment[p][d] = Auxils.adjustToRange(environment[p][d], comm.minEnv, comm.maxEnv);
                    adjustFitness(p, d);
                }
            }
        }
    } else {
        for (int d = 0; d < comm.envDims; d++) {
            for (int p = 0; p < comm.nbrPatches; p++) {
                if (Auxils.random.nextDouble() <= comm.pChange[pc]) {
                    step = Auxils.gaussianSampler.sample() * comm.envStep[esPos];
                    environment[p][d] += step;
                    environment[p][d] = Auxils.adjustToRange(environment[p][d], comm.minEnv, comm.maxEnv);
                    adjustFitness(p, d);
                }
            }
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

    void mortality() {
        double fit = 0.;
        // double surv = 0.;
        int p;

        for (int i = 0; i < totSites; i++) {
            p = patch[i];
            if (alive[i]) {
// soft selection
                fit = (fitness[i] / maxFitness[p]);
// hard selection
                // fit = fitness[i];
// regular selection
                // alive[i] = Auxils.random.nextDouble() < (1 - comm.d) * fit;

                if (Auxils.random.nextDouble() >= (1 - comm.d) * fit)
                    removeInd(i);
// density dependent selection
                // surv = Math.max(1 + popS[p]/(1*comm.microsites)*(fit-1), 0);
                // alive[i] = Auxils.random.nextDouble() < (1 - comm.d) * surv;
// K dependent selection
                // surv = Math.max(1 - popS[p]/(fit*2*comm.microsites), 0);
                // alive[i] = Auxils.random.nextDouble() < (1 - comm.d) * surv;
// dens2 selection
                // surv = Math.max(fit*(1 - popS[p]/(comm.microsites*2)), 0);
                // alive[i] = Auxils.random.nextDouble() < (1 - comm.d) * surv;
            }
        }
    }
    
    // dispersal
    void disperse() {
        nbrDisp = 0;
        // int aliveDisp = 0;
        double maxDisp = 0;
        double iDisp;
        int maxN = 0;
        double[] pEmpty = new double[comm.nbrPatches];
        maxDisp = Auxils.arrayMax(pDispSum);
        maxN = Auxils.arrayMax(popN);
        
        for (int p = 0; p < comm.nbrPatches; p++) {
            if (evol.mutationRateDisp == 0)
                pEmpty[p] = (1. + Math.round(maxN*comm.dispRate[drPos] - popN[p]*comm.dispRate[drPos]))/(double)(comm.microsites - popN[p]);
            else
                pEmpty[p] = (1. + Math.round(maxDisp - pDispSum[p]))/(double)(comm.microsites - popN[p]);
        }

        for (int i = 0; i < totSites; i++) {
            // if(alive[i] && Auxils.random.nextDouble() < comm.dispRate[drPos]) {
            //     posDisp[nbrDisp++] = i;
            // }
            if (alive[i]) {
                if (evol.mutationRateDisp == 0)
                    iDisp = comm.dispRate[drPos];
                else
                    iDisp = pDisp[i];
                if (Auxils.random.nextDouble() < iDisp)
                    posDisp[nbrDisp++] = i;
            } else if (Auxils.random.nextDouble() < pEmpty[patch[i]])
                posDisp[nbrDisp++] = i;
        }

        if(nbrDisp > 0) {
            // if (nbrDisp > aliveDisp) {
            //     System.out.println("     disp: " + nbrDisp + ",  alive disp: " + aliveDisp + ",  popsize: " + metapopSize());
            // }
            int oldPos, newPos, i2;
            // int[] dispShuffle = Arrays.copyOf(posDisp, nbrDisp);
            int[] dispShuffle = new int[nbrDisp];
            System.arraycopy(posDisp, 0, dispShuffle, 0, nbrDisp);
            Auxils.arrayShuffle(dispShuffle);
            // byte[] tempGen = Arrays.copyOf(genotype[dispShuffle[0]], 2 * evol.allLoci);
            byte[] tempGen = new byte[2 * evol.allLoci];
            System.arraycopy(genotype[dispShuffle[0]], 0, tempGen, 0, 2 * evol.allLoci);
            boolean tempAlive = alive[dispShuffle[0]];
            if (tempAlive)
                removeInd(dispShuffle[0]);
            for (int i = 1; i < nbrDisp; i++) {
                oldPos = dispShuffle[i];
                newPos = dispShuffle[i - 1];
                i2 = i+1;
                while(patch[oldPos] == patch[newPos] && i2 < nbrDisp) {
                    dispShuffle[i] = dispShuffle[i2];
                    dispShuffle[i2] = oldPos;
                    oldPos = dispShuffle[i];
                    i2++;
                }
                // System.out.println("     disp: ");
                // System.out.println("     old patch = " + patch[oldPos] + "; old pos = " + oldPos);
                // System.out.println("     new patch = " + patch[newPos] + "; new pos = " + newPos);
                if (alive[oldPos]) {
                    System.arraycopy(genotype[oldPos], 0, genotype[newPos], 0, 2 * evol.allLoci);
                    settleRest(newPos, oldPos);
                    removeInd(oldPos);
                } else {
                    alive[newPos] = false;
                }
            }
            newPos = dispShuffle[nbrDisp - 1];
            if (tempAlive) {
                System.arraycopy(tempGen, 0, genotype[newPos], 0, 2 * evol.allLoci);
                settleRest(newPos, dispShuffle[0]);
            } else {
                alive[newPos] = false;
            }
        }
    }

    void contributionAdults() {
        Arrays.fill(endPosFathers, 0);
        Arrays.fill(endPosMothers, 0);
        Arrays.fill(nbrEmpty, 0);
        Arrays.fill(production, 0.);

        double contr = 0.;
        int p;

        // double fit;

        // int[] popS = new int[comm.nbrPatches];
        // Arrays.fill(popS, 0);
        // for (p = 0; p < comm.nbrPatches; p++) {
        //     popS[p] = popSize(p);
        // }

        for (int i = 0; i < totSites; i++) {
            p = patch[i];
            if (alive[i]) {

//                for (int l = 0; l < (2*evol.allLoci); l++) {
//                    migrationGenotype[i][l] += 1;
//                }
                
                // fit = (fitness[i] / maxFitness[p]);

                contr = 1;
                sexAdults[i] = Auxils.random.nextDouble() <= pSex[i];
                // sexAdults[i] = Auxils.random.nextDouble() <= (1 - fit)*pSex[i];
                if (sexAdults[i]) {
                    // fathersPos[p][endPosFathers[p]] = i;
                    // fathersProb[p][endPosFathers[p]++] = contr;
                    fathersPos[p][endPosFathers[p]++] = i;
                    contr *= comm.demogrCost[dcPos];
                }
                mothersPos[p][endPosMothers[p]] = i;
                mothersProb[p][endPosMothers[p]++] = contr;
            } else {
                contr = 0.;
                posEmpty[p][nbrEmpty[p]++] = i;
            }
        }
        for (p = 0; p < comm.nbrPatches; p++) {
            if (endPosMothers[p] > 0) {
                // mothersCumProb[p] = Arrays.copyOf(mothersProb[p], endPosMothers[p]);

                // mothersCumProb[p] = new double[endPosMothers[p]];
                System.arraycopy(mothersProb[p], 0, mothersCumProb[p], 0, endPosMothers[p]);

                Auxils.arrayCumSum(mothersCumProb[p], endPosMothers[p]);
                production[p] = mothersCumProb[p][endPosMothers[p] - 1] * comm.r;
                Auxils.arrayDiv(mothersCumProb[p], endPosMothers[p], mothersCumProb[p][endPosMothers[p] - 1]);
                // if (endPosFathers[p] > 0) {
                //     // fathersCumProb[p] = Arrays.copyOf(fathersProb[p], endPosFathers[p]);

                //     // fathersCumProb[p] = new double[endPosFathers[p]];
                //     System.arraycopy(fathersProb[p], 0, fathersCumProb[p], 0, endPosFathers[p]);
    
                //         Auxils.arrayCumSum(fathersCumProb[p], endPosFathers[p]);
                //     Auxils.arrayDiv(fathersCumProb[p], endPosFathers[p], fathersCumProb[p][endPosFathers[p] - 1]);
                // }
            }
        }
//        System.out.println("      motherscumprob = ");
//        System.out.println(Arrays.deepToString(mothersCumProb));
    }

    void reproduction() {
        // int[] posOffspring;
        int m, f;
        double prod;

        for (int p = 0; p < comm.nbrPatches; p++) {
            if (production[p] > 0.) {
                prod = (production[p] < 1) ? ((Auxils.random.nextDouble() < production[p]) ? 1. : 0.) : production[p];
                nbrSettled = Math.min(nbrEmpty[p], (int) prod);
                // posOffspring = Auxils.arraySample(nbrSettled, Arrays.copyOf(posEmpty[p], nbrEmpty[p]));
                int[] posOffspring = Auxils.arraySample(nbrSettled, posEmpty[p], nbrEmpty[p]);
                //sampling parents with replacement!
                for (int i = 0; i < nbrSettled; i++) {
                    m = mothersPos[p][Auxils.randIntCumProb(mothersCumProb[p], endPosMothers[p])];
                    if (sexAdults[m]) {
                        //selfing allowed!
                        // f = fathersPos[p][Auxils.randIntCumProb(fathersCumProb[p], endPosFathers[p])];
                        f = fathersPos[p][Auxils.random.nextInt(endPosFathers[p])];

                        while (pSex[f] != pSex[m]) {
                            // f = fathersPos[p][Auxils.randIntCumProb(fathersCumProb[p], endPosFathers[p])];
                            f = fathersPos[p][Auxils.random.nextInt(endPosFathers[p])];
                        }

                        settle(posOffspring[i], m, f);
                        //selfing not allowed!
                        // if(endPosFathers[patchMother] > 1) {
                        //     f = fathersPos[patchMother][Auxils.randIntCumProb(fathersCumProb[patchMother])];
                        //     while (f == m) {
                        //         f = fathersPos[patchMother][Auxils.randIntCumProb(fathersCumProb[patchMother])];
                        //     }
                        //     settle(posOffspring[i], m, f);
                        // }
                    } else {
                        settle(posOffspring[i], m);
                    }
                }
            }
        }
    }

    /* install newborns and inherit traits from the parent(s)
     * including mutation */
    void settle(int pos, int m, int f) {
        inherit(pos, m, f);
        mutate(pos);
        settleRest(pos, m);
    }

    void settle(int pos, int m) {
        inherit(pos, m);
        mutate(pos);
        settleRest(pos, m);
    }

    void settleRest(int pos, int m) {
        int p = patch[pos];
        alive[pos] = true;
        popN[p]++;
        if (patch[m] != p) {
            newMigrant(pos);
//            for (int l = 0; l < evol.allLoci; l++) {
////                migrationGenotype[pos][evol.allMother[l]] = 1;
////                migrationGenotype[pos][evol.allFather[l]] = 1;
//                migrationGenotype[pos][evol.allMother[l]] = migrationCounter;
//                migrationGenotype[pos][evol.allFather[l]] = migrationCounter;
//            }
//            migrationCounter++;
        }
        fitness[pos] = 1;
        for (int tr = 0; tr < comm.traits; tr++) {
            traitPhenotype[pos][tr] = calcPhenotype(pos, tr);
            traitFitness[pos][tr] = calcFitness(traitPhenotype[pos][tr], environment[p][comm.traitDim[tr]]);
            fitness[pos] *= traitFitness[pos][tr];
        }
        // if (maxFitness[p] < fitness[pos])
        //     maxFitness[p] = fitness[pos];
        // pSex[pos] = Math.min(1, Math.max(0, Auxils.arrayMean(Auxils.arrayElements(genotype[pos], evol.sexGenes))));
        pSex[pos] = Math.min(1, Math.max(0, Auxils.arrayMean(genotype[pos], evol.sexGenes)));
//        pDisp[pos] = Math.pow(Math.min(1, Math.max(0, Auxils.arrayMean(Auxils.arrayElements(genotype[pos], evol.dispGenes)))), evol.dispPow) / evol.dispCorr;
        // pDisp[pos] = Math.exp(Math.log(evol.refDisp)*(1 - Math.min(1, Math.max(0, Auxils.arrayMean(Auxils.arrayElements(genotype[pos], evol.dispGenes)))))/(1 - 0.2));
        pDisp[pos] = Math.exp(Math.log(evol.refDisp)*(1 - Math.min(1, Math.max(0, Auxils.arrayMean(genotype[pos], evol.dispGenes))))/(1 - 0.2));

        pDispSum[p] += pDisp[pos];
    }

    void removeInd (int pos) {
        int p = patch[pos];
        alive[pos] = false;
        popN[p]--;
        pDispSum[p] -= pDisp[pos];
    }


//    void settle(int p, int[] posOffspring, int[] patchOrigin) {
//        int pos;
//        for (int i = 0; i < nbrSettled; i++) {
//            pos = posOffspring[i];
//            alive[pos] = true;
//            System.arraycopy(newborns[p][i], 0, genotype[pos], 0, 2 * evol.allLoci);
//
//            if (patchOrigin[i] != p) {
//                for (int l = 0; l < (2*evol.allLoci); l++) {
//                    migrationGenotype[pos][l] = 1;
//                }
//            } else
//                System.arraycopy(migrationNewborns[p][i], 0, migrationGenotype[pos], 0, 2 * evol.allLoci);
//
//
//            fitness[pos] = 1;
//            for (int tr = 0; tr < comm.traits; tr++) {
//                traitPhenotype[pos][tr] = calcPhenotype(pos, tr);
//                traitFitness[pos][tr] = calcFitness(traitPhenotype[pos][tr], environment[p][comm.traitDim[tr]]);
//                fitness[pos] *= traitFitness[pos][tr];
//            }
//            if (maxFitness[p] < fitness[pos])
//                maxFitness[p] = fitness[pos];
//            pSex[pos] = Math.min(1, Math.max(0, Auxils.arrayMean(Auxils.arrayElements(genotype[pos], evol.sexGenes))));
//        }
//    }

    /* inheritance for asexual reproduction (one parent) */
    void inherit(int posOffspring, int posParent) {
        System.arraycopy(genotype[posParent], 0, genotype[posOffspring], 0, 2 * evol.allLoci);

        System.arraycopy(migrationGenotype[posParent], 0, migrationGenotype[posOffspring], 0, 2 * evol.allLoci);
//        for (int l = 0; l < (2*evol.allLoci); l++) {
//            migrationNewborns[p][posOffspring][l] += 1;
//        }
    }

//    void inherit(int p, int posOffspring, int posParent) {
//        System.arraycopy(genotype[posParent], 0, newborns[p][posOffspring], 0, 2 * evol.allLoci);
//
//        System.arraycopy(migrationGenotype[posParent], 0, migrationNewborns[p][posOffspring], 0, 2 * evol.allLoci);
////        for (int l = 0; l < (2*evol.allLoci); l++) {
////            migrationNewborns[p][posOffspring][l] += 1;
////        }
//    }

    /* inheritance for sexual reproduction (two parent) */
    void inherit_1(int posOffspring, int posMother, int posFather) {
        for (int l = 0; l < evol.allLoci; l++) {
//            newborns[p][posOffspring][evol.allMother[l]] = genotype[posMother][Auxils.random.nextBoolean() ? evol.allMother[l] : evol.allFather[l]];
//            newborns[p][posOffspring][evol.allFather[l]] = genotype[posFather][Auxils.random.nextBoolean() ? evol.allMother[l] : evol.allFather[l]];

            if(Auxils.random.nextBoolean()) {
                genotype[posOffspring][evol.allMother[l]] = genotype[posMother][evol.allMother[l]];
//                migrationNewborns[p][posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allMother[l]] + 1;
                migrationGenotype[posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allMother[l]];
            } else {
                genotype[posOffspring][evol.allMother[l]] = genotype[posMother][evol.allFather[l]];
//                migrationNewborns[p][posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allFather[l]] + 1;
                migrationGenotype[posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allFather[l]];
            }
            if(Auxils.random.nextBoolean()) {
                genotype[posOffspring][evol.allFather[l]] = genotype[posFather][evol.allMother[l]];
//                migrationNewborns[p][posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allMother[l]] + 1;
                migrationGenotype[posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allMother[l]];
            } else {
                genotype[posOffspring][evol.allFather[l]] = genotype[posFather][evol.allFather[l]];
//                migrationNewborns[p][posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allFather[l]] + 1;
                migrationGenotype[posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allFather[l]];
            }
        }
    }

    /* inheritance for sexual reproduction (two parent) */
    void inherit_2(int posOffspring, int posMother, int posFather) {
        int kmf, kff;
        int[] mfLocs, ffLocs;
        CombinationSampler mfCombinationSampler, ffCombinationSampler;

        kmf = Auxils.binomialSamplerInherit.sample();
        mfCombinationSampler = new CombinationSampler(Auxils.random, evol.allLoci, kmf);
        mfLocs = mfCombinationSampler.sample();

        System.arraycopy(genotype[posMother], 0, genotype[posOffspring], 0, evol.allLoci);
        System.arraycopy(migrationGenotype[posMother], 0, migrationGenotype[posOffspring], 0, evol.allLoci);
        for (int l : mfLocs) {
            genotype[posOffspring][evol.allMother[l]] = genotype[posMother][evol.allFather[l]];
            migrationGenotype[posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allFather[l]];
        }

        kff = Auxils.binomialSamplerInherit.sample();
        ffCombinationSampler = new CombinationSampler(Auxils.random, evol.allLoci, kff);
        ffLocs = ffCombinationSampler.sample();
        System.arraycopy(genotype[posFather], 0, genotype[posOffspring], evol.allLoci, evol.allLoci);
        System.arraycopy(migrationGenotype[posFather], 0, migrationGenotype[posOffspring], evol.allLoci, evol.allLoci);
        for (int l : ffLocs) {
            genotype[posOffspring][evol.allFather[l]] = genotype[posFather][evol.allFather[l]];
            migrationGenotype[posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allFather[l]];
        }
    }

    void inherit_3(int posOffspring, int posMother, int posFather) {
        int kmf, kff, l;
        int[] shuffleLocs = Auxils.enumArray(0, evol.allLoci - 1);

        Auxils.arrayShuffle(shuffleLocs);
        kmf = Auxils.binomialSamplerInherit.sample();
        for (int i = 0; i < kmf; i++) {
            l = shuffleLocs[i];
            genotype[posOffspring][evol.allMother[l]] = genotype[posMother][evol.allMother[l]];
            migrationGenotype[posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allMother[l]];
        }
        for (int i = kmf; i < evol.allLoci; i++) {
            l = shuffleLocs[i];
            genotype[posOffspring][evol.allMother[l]] = genotype[posMother][evol.allFather[l]];
            migrationGenotype[posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allFather[l]];
        }

        Auxils.arrayShuffle(shuffleLocs);
        kff = Auxils.binomialSamplerInherit.sample();
        for (int i = 0; i < kff; i++) {
            l = shuffleLocs[i];
            genotype[posOffspring][evol.allFather[l]] = genotype[posFather][evol.allMother[l]];
            migrationGenotype[posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allMother[l]];
        }
        for (int i = kff; i < evol.allLoci; i++) {
            l = shuffleLocs[i];
            genotype[posOffspring][evol.allFather[l]] = genotype[posFather][evol.allFather[l]];
            migrationGenotype[posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allFather[l]];
        }
    }

    void inherit_4(int posOffspring, int posMother, int posFather) {
        int kmf, kff, l;

        Auxils.arrayShuffle(evol.shuffleLocs);
        kmf = Auxils.binomialSamplerInherit.sample();
        for (int i = 0; i < kmf; i++) {
            l = evol.shuffleLocs[i];
            genotype[posOffspring][evol.allMother[l]] = genotype[posMother][evol.allMother[l]];
            migrationGenotype[posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allMother[l]];
        }
        for (int i = kmf; i < evol.allLoci; i++) {
            l = evol.shuffleLocs[i];
            genotype[posOffspring][evol.allMother[l]] = genotype[posMother][evol.allFather[l]];
            migrationGenotype[posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allFather[l]];
        }

        Auxils.arrayShuffle(evol.shuffleLocs);
        kff = Auxils.binomialSamplerInherit.sample();
        for (int i = 0; i < kff; i++) {
            l = evol.shuffleLocs[i];
            genotype[posOffspring][evol.allFather[l]] = genotype[posFather][evol.allMother[l]];
            migrationGenotype[posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allMother[l]];
        }
        for (int i = kff; i < evol.allLoci; i++) {
            l = evol.shuffleLocs[i];
            genotype[posOffspring][evol.allFather[l]] = genotype[posFather][evol.allFather[l]];
            migrationGenotype[posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allFather[l]];
        }
    }

    void inherit_5(int posOffspring, int posMother, int posFather) {
        int kmf, kff, l, lm, lf;

        kmf = Auxils.binomialSamplerInherit.sample();
        Auxils.arrayShuffle(evol.shuffleLocs, kmf);
        System.arraycopy(genotype[posMother], 0, genotype[posOffspring], 0, evol.allLoci);
        System.arraycopy(migrationGenotype[posMother], 0, migrationGenotype[posOffspring], 0, evol.allLoci);
        // for (int i = 0; i < kmf; i++) {
        //     l = evol.shuffleLocs[i];
        //     genotype[posOffspring][evol.allMother[l]] = genotype[posMother][evol.allMother[l]];
        //     migrationGenotype[posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allMother[l]];
        // }
        for (int i = 0; i < kmf; i++) {
            l = evol.shuffleLocs[i];
            lm = evol.allMother[l];
            lf = evol.allFather[l];
            genotype[posOffspring][lm] = genotype[posMother][lf];
            migrationGenotype[posOffspring][lm] = migrationGenotype[posMother][lf];
        }

        kff = Auxils.binomialSamplerInherit.sample();
        Auxils.arrayShuffle(evol.shuffleLocs, kff);
        System.arraycopy(genotype[posFather], 0, genotype[posOffspring], evol.allLoci, evol.allLoci);
        System.arraycopy(migrationGenotype[posFather], 0, migrationGenotype[posOffspring], evol.allLoci, evol.allLoci);
        // for (int i = 0; i < kff; i++) {
        //     l = evol.shuffleLocs[i];
        //     genotype[posOffspring][evol.allFather[l]] = genotype[posFather][evol.allMother[l]];
        //     migrationGenotype[posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allMother[l]];
        // }
        for (int i = 0; i < kff; i++) {
            l = evol.shuffleLocs[i];
            lf = evol.allFather[l];
            genotype[posOffspring][lf] = genotype[posFather][lf];
            migrationGenotype[posOffspring][lf] = migrationGenotype[posFather][lf];
        }
    }

    void inherit(int posOffspring, int posMother, int posFather) {
        int kmf, kff, l, lm, lf, index;

        kmf = Auxils.binomialSamplerInherit.sample();
        // for (int i = 0; i < kmf; i++) {
        //     index = Auxils.random.nextInt(i, evol.allLoci);
        //     l = evol.shuffleLocs[index];
        //     evol.shuffleLocs[index] = evol.shuffleLocs[i];
        //     evol.shuffleLocs[i] = l;
        //     lm = evol.allMother[l];
        //     // lf = evol.allFather[l];
        //     genotype[posOffspring][lm] = genotype[posMother][lm];
        //     migrationGenotype[posOffspring][lm] = migrationGenotype[posMother][lm];
        // }
        // for (int i = kmf; i < evol.allLoci; i++) {
        //     l = evol.shuffleLocs[i];
        //     lm = evol.allMother[l];
        //     lf = evol.allFather[l];
        //     genotype[posOffspring][lm] = genotype[posMother][lf];
        //     migrationGenotype[posOffspring][lm] = migrationGenotype[posMother][lf];
        // }
        System.arraycopy(genotype[posMother], 0, genotype[posOffspring], 0, evol.allLoci);
        System.arraycopy(migrationGenotype[posMother], 0, migrationGenotype[posOffspring], 0, evol.allLoci);
        for (int i = 0; i < kmf; i++) {
            index = Auxils.random.nextInt(i, evol.allLoci);
            l = evol.shuffleLocs[index];
            evol.shuffleLocs[index] = evol.shuffleLocs[i];
            evol.shuffleLocs[i] = l;
            lm = evol.allMother[l];
            lf = evol.allFather[l];
            genotype[posOffspring][lm] = genotype[posMother][lf];
            migrationGenotype[posOffspring][lm] = migrationGenotype[posMother][lf];
        }

        kff = Auxils.binomialSamplerInherit.sample();
        // for (int i = 0; i < kff; i++) {
        //     index = Auxils.random.nextInt(i, evol.allLoci);
        //     l = evol.shuffleLocs[index];
        //     evol.shuffleLocs[index] = evol.shuffleLocs[i];
        //     evol.shuffleLocs[i] = l;
        //     lm = evol.allMother[l];
        //     lf = evol.allFather[l];
        //     genotype[posOffspring][lf] = genotype[posFather][lm];
        //     migrationGenotype[posOffspring][lf] = migrationGenotype[posFather][lm];
        // }
        // for (int i = kff; i < evol.allLoci; i++) {
        //     l = evol.shuffleLocs[i];
        //     lf = evol.allFather[l];
        //     genotype[posOffspring][lf] = genotype[posFather][lf];
        //     migrationGenotype[posOffspring][lf] = migrationGenotype[posFather][lf];
        // }
        System.arraycopy(genotype[posFather], 0, genotype[posOffspring], evol.allLoci, evol.allLoci);
        System.arraycopy(migrationGenotype[posFather], 0, migrationGenotype[posOffspring], evol.allLoci, evol.allLoci);
        for (int i = 0; i < kff; i++) {
            index = Auxils.random.nextInt(i, evol.allLoci);
            l = evol.shuffleLocs[index];
            evol.shuffleLocs[index] = evol.shuffleLocs[i];
            evol.shuffleLocs[i] = l;
            lf = evol.allFather[l];
            genotype[posOffspring][lf] = genotype[posFather][lf];
            migrationGenotype[posOffspring][lf] = migrationGenotype[posFather][lf];
        }
    }

    void inherit_7(int posOffspring, int posMother, int posFather) {
        int kmf, kff, l, lm, lf, index;

        kmf = Auxils.binomialSamplerInherit.sample();
        // for (int i = 0; i < kmf; i++) {
        //     index = Auxils.random.nextInt(i, evol.allLoci);
        //     l = evol.shuffleLocs[index];
        //     evol.shuffleLocs[index] = evol.shuffleLocs[i];
        //     evol.shuffleLocs[i] = l;
        //     lm = evol.allMother[l];
        //     // lf = evol.allFather[l];
        //     genotype[posOffspring][lm] = genotype[posMother][lm];
        //     migrationGenotype[posOffspring][lm] = migrationGenotype[posMother][lm];
        // }
        // for (int i = kmf; i < evol.allLoci; i++) {
        //     l = evol.shuffleLocs[i];
        //     lm = evol.allMother[l];
        //     lf = evol.allFather[l];
        //     genotype[posOffspring][lm] = genotype[posMother][lf];
        //     migrationGenotype[posOffspring][lm] = migrationGenotype[posMother][lf];
        // }
        System.arraycopy(genotype[posMother], 0, genotype[posOffspring], 0, evol.allLoci);
        System.arraycopy(migrationGenotype[posMother], 0, migrationGenotype[posOffspring], 0, evol.allLoci);
        for (int i = evol.allLoci - 1; i > evol.allLoci - kmf - 1; i--) {
            index = Auxils.random.nextInt(i + 1);
            l = evol.shuffleLocs[index];
            evol.shuffleLocs[index] = evol.shuffleLocs[i];
            evol.shuffleLocs[i] = l;
            lm = evol.allMother[l];
            lf = evol.allFather[l];
            genotype[posOffspring][lm] = genotype[posMother][lf];
            migrationGenotype[posOffspring][lm] = migrationGenotype[posMother][lf];
        }

        kff = Auxils.binomialSamplerInherit.sample();
        // for (int i = 0; i < kff; i++) {
        //     index = Auxils.random.nextInt(i, evol.allLoci);
        //     l = evol.shuffleLocs[index];
        //     evol.shuffleLocs[index] = evol.shuffleLocs[i];
        //     evol.shuffleLocs[i] = l;
        //     lm = evol.allMother[l];
        //     lf = evol.allFather[l];
        //     genotype[posOffspring][lf] = genotype[posFather][lm];
        //     migrationGenotype[posOffspring][lf] = migrationGenotype[posFather][lm];
        // }
        // for (int i = kff; i < evol.allLoci; i++) {
        //     l = evol.shuffleLocs[i];
        //     lf = evol.allFather[l];
        //     genotype[posOffspring][lf] = genotype[posFather][lf];
        //     migrationGenotype[posOffspring][lf] = migrationGenotype[posFather][lf];
        // }
        System.arraycopy(genotype[posFather], 0, genotype[posOffspring], evol.allLoci, evol.allLoci);
        System.arraycopy(migrationGenotype[posFather], 0, migrationGenotype[posOffspring], evol.allLoci, evol.allLoci);
        for (int i = evol.allLoci - 1; i > evol.allLoci - kff - 1; i--) {
            index = Auxils.random.nextInt(i + 1);
            l = evol.shuffleLocs[index];
            evol.shuffleLocs[index] = evol.shuffleLocs[i];
            evol.shuffleLocs[i] = l;
            lf = evol.allFather[l];
            genotype[posOffspring][lf] = genotype[posFather][lf];
            migrationGenotype[posOffspring][lf] = migrationGenotype[posFather][lf];
        }
    }

    void mutate_1(int posOffspring) {
        int k, l;
        int tempAll;
        // int[] somMutLocs;

        k = Auxils.binomialSamplerSomatic.sample();
        if (k > 0) {
            CombinationSampler combinationSampler = new CombinationSampler(Auxils.random, evol.traitLoci*2, k);
            // somMutLocs = Auxils.arrayElements(evol.somGenes, combinationSampler.sample());
            int[] somMutLocs = combinationSampler.sample();
            for (int i : somMutLocs) {
                l = evol.somGenes[i];
                genotype[posOffspring][l] += (Auxils.random.nextBoolean() ? -1 : 1);
            }
        }

        if (comm.sexType.equals("SWITCH")) {
            if (Auxils.random.nextDouble() <= evol.mutationRateSex) {
                for (int i : evol.sexGenes) {
                    genotype[posOffspring][i] = (byte) ((genotype[posOffspring][i] == 0) ? 1 : 0);
                }
            }
        } else {
            // int[] sexMutLocs;
            double pSexTemp;
            k = Auxils.binomialSamplerSex.sample();
            if (k > 0) {
                pSexTemp = Auxils.arrayMean(genotype[posOffspring], evol.sexGenes);
                CombinationSampler combinationSampler = new CombinationSampler(Auxils.random, evol.sexLoci * 2, k);
                // sexMutLocs = Auxils.arrayElements(evol.sexGenes, combinationSampler.sample());
                int[] sexMutLocs = combinationSampler.sample();
                for (int i : sexMutLocs) {
                    l = evol.sexGenes[i];
                    if (pSexTemp <= 0.) {
                        genotype[posOffspring][l] += 1;
                        // pSexTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.sexGenes));
                        pSexTemp += 1./evol.sexGenes.length;
                    }
                    else if (pSexTemp >= 1.) {
                        genotype[posOffspring][l] -= 1;
                        // pSexTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.sexGenes));
                        pSexTemp -= 1./evol.sexGenes.length;
                    }
                    else {
                        tempAll = (Auxils.random.nextBoolean() ? -1 : 1);
                        genotype[posOffspring][l] += tempAll;
                        // pSexTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.sexGenes));
                        pSexTemp += ((double)tempAll)/evol.sexGenes.length;
                    }
                }
            }
        }

        // int[] dispMutLocs;
        double pDispTemp;
        k = Auxils.binomialSamplerDisp.sample();
        if (k > 0) {
            pDispTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.dispGenes));
            CombinationSampler combinationSampler = new CombinationSampler(Auxils.random, evol.dispLoci * 2, k);
            // int[] dispMutLocs = Auxils.arrayElements(evol.dispGenes, combinationSampler.sample());
            int[] dispMutLocs = combinationSampler.sample();
            for (int i : dispMutLocs) {
                l = evol.dispGenes[i];
                if (pDispTemp <= 0.) {
                    genotype[posOffspring][l] += 1;
                    // pDispTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.dispGenes));
                    pDispTemp += 1./evol.dispGenes.length;
                }
                else if (pDispTemp >= 1.) {
                    genotype[posOffspring][l] -= 1;
                    // pDispTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.dispGenes));
                    pDispTemp -= 1./evol.dispGenes.length;
                }
                else {
                    tempAll = (Auxils.random.nextBoolean() ? -1 : 1);
                    // genotype[posOffspring][l] += (Auxils.random.nextBoolean() ? -1 : 1);
                    genotype[posOffspring][l] += tempAll;
                    // pDispTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.dispGenes));
                    pDispTemp += ((double)tempAll)/evol.dispGenes.length;
                }
            }
        }
    }

    void mutate_2(int posOffspring) {
        int k, l;
        int tempAll;
        // int[] somMutLocs;

        k = Auxils.binomialSamplerSomatic.sample();
        if (k > 0) {
            Auxils.arrayShuffle(evol.shuffleTrait, k);
            for (int i = 0; i < k; i++) {
                l = evol.somGenes[evol.shuffleTrait[i]];
                genotype[posOffspring][l] += (Auxils.random.nextBoolean() ? -1 : 1);
            }
        }

        if (comm.sexType.equals("SWITCH")) {
            if (Auxils.random.nextDouble() <= evol.mutationRateSex) {
                for (int i : evol.sexGenes) {
                    genotype[posOffspring][i] = (byte) ((genotype[posOffspring][i] == 0) ? 1 : 0);
                }
            }
        } else {
            // int[] sexMutLocs;
            double pSexTemp;
            k = Auxils.binomialSamplerSex.sample();
            if (k > 0) {
                pSexTemp = Auxils.arrayMean(genotype[posOffspring], evol.sexGenes);
                Auxils.arrayShuffle(evol.shuffleSex, k);
                for (int i = 0; i < k; i++) {
                    l = evol.sexGenes[evol.shuffleSex[i]];
                    if (pSexTemp <= 0.) {
                        genotype[posOffspring][l] += 1;
                        // pSexTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.sexGenes));
                        pSexTemp += 1./evol.sexGenes.length;
                    }
                    else if (pSexTemp >= 1.) {
                        genotype[posOffspring][l] -= 1;
                        // pSexTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.sexGenes));
                        pSexTemp -= 1./evol.sexGenes.length;
                    }
                    else {
                        tempAll = (Auxils.random.nextBoolean() ? -1 : 1);
                        genotype[posOffspring][l] += tempAll;
                        // pSexTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.sexGenes));
                        pSexTemp += ((double)tempAll)/evol.sexGenes.length;
                    }
                }
            }
        }

        // int[] dispMutLocs;
        double pDispTemp;
        k = Auxils.binomialSamplerDisp.sample();
        if (k > 0) {
            pDispTemp = Auxils.arrayMean(genotype[posOffspring], evol.dispGenes);
            Auxils.arrayShuffle(evol.shuffleDisp, k);
            for (int i = 0; i < k; i++) {
                l = evol.dispGenes[evol.shuffleDisp[i]];
                if (pDispTemp <= 0.) {
                    genotype[posOffspring][l] += 1;
                    // pDispTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.dispGenes));
                    pDispTemp += 1./evol.dispGenes.length;
                }
                else if (pDispTemp >= 1.) {
                    genotype[posOffspring][l] -= 1;
                    // pDispTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.dispGenes));
                    pDispTemp -= 1./evol.dispGenes.length;
                }
                else {
                    tempAll = (Auxils.random.nextBoolean() ? -1 : 1);
                    // genotype[posOffspring][l] += (Auxils.random.nextBoolean() ? -1 : 1);
                    genotype[posOffspring][l] += tempAll;
                    // pDispTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.dispGenes));
                    pDispTemp += ((double)tempAll)/evol.dispGenes.length;
                }
            }
        }
    }

    void mutate(int posOffspring) {
        int k, l, index;
        int tempAll;
        // int[] somMutLocs;

        k = Auxils.binomialSamplerSomatic.sample();
        if (k > 0) {
            // Auxils.arrayShuffle(evol.shuffleTrait, k);
            for (int i = 0; i < k; i++) {
                index = Auxils.random.nextInt(i, evol.shuffleTrait.length);
                l = evol.shuffleTrait[index];
                evol.shuffleTrait[index] = evol.shuffleTrait[i];
                evol.shuffleTrait[i] = l;
                genotype[posOffspring][evol.somGenes[l]] += (Auxils.random.nextBoolean() ? -1 : 1);
            }
        }

        if (comm.sexType.equals("SWITCH")) {
            if (Auxils.random.nextDouble() <= evol.mutationRateSex) {
                for (int i : evol.sexGenes) {
                    genotype[posOffspring][i] = (byte) ((genotype[posOffspring][i] == 0) ? 1 : 0);
                }
            }
        } else {
            // int[] sexMutLocs;
            double pSexTemp;
            k = Auxils.binomialSamplerSex.sample();
            if (k > 0) {
                pSexTemp = Auxils.arrayMean(genotype[posOffspring], evol.sexGenes);
                // Auxils.arrayShuffle(evol.shuffleSex, k);
                for (int i = 0; i < k; i++) {
                    index = Auxils.random.nextInt(i, evol.shuffleSex.length);
                    l = evol.shuffleSex[index];
                    evol.shuffleSex[index] = evol.shuffleSex[i];
                    evol.shuffleSex[i] = l;
                    if (pSexTemp <= 0.) {
                        genotype[posOffspring][evol.sexGenes[l]] += 1;
                        // pSexTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.sexGenes));
                        pSexTemp += 1./evol.sexGenes.length;
                    }
                    else if (pSexTemp >= 1.) {
                        genotype[posOffspring][evol.sexGenes[l]] -= 1;
                        // pSexTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.sexGenes));
                        pSexTemp -= 1./evol.sexGenes.length;
                    }
                    else {
                        tempAll = (Auxils.random.nextBoolean() ? -1 : 1);
                        genotype[posOffspring][evol.sexGenes[l]] += tempAll;
                        // pSexTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.sexGenes));
                        pSexTemp += ((double)tempAll)/evol.sexGenes.length;
                    }
                }
            }
        }

        // int[] dispMutLocs;
        double pDispTemp;
        k = Auxils.binomialSamplerDisp.sample();
        if (k > 0) {
            pDispTemp = Auxils.arrayMean(genotype[posOffspring], evol.dispGenes);
            // Auxils.arrayShuffle(evol.shuffleDisp, k);
            for (int i = 0; i < k; i++) {
                index = Auxils.random.nextInt(i, evol.shuffleDisp.length);
                l = evol.shuffleDisp[index];
                evol.shuffleDisp[index] = evol.shuffleDisp[i];
                evol.shuffleDisp[i] = l;
                if (pDispTemp <= 0.) {
                    genotype[posOffspring][evol.dispGenes[l]] += 1;
                    // pDispTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.dispGenes));
                    pDispTemp += 1./evol.dispGenes.length;
                }
                else if (pDispTemp >= 1.) {
                    genotype[posOffspring][evol.dispGenes[l]] -= 1;
                    // pDispTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.dispGenes));
                    pDispTemp -= 1./evol.dispGenes.length;
                }
                else {
                    tempAll = (Auxils.random.nextBoolean() ? -1 : 1);
                    // genotype[posOffspring][l] += (Auxils.random.nextBoolean() ? -1 : 1);
                    genotype[posOffspring][evol.dispGenes[l]] += tempAll;
                    // pDispTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring], evol.dispGenes));
                    pDispTemp += ((double)tempAll)/evol.dispGenes.length;
                }
            }
        }
    }

    void newMigrant(int pos) {
        int p = patch[pos];
        Arrays.fill(migrationGenotype[pos], migrationCounter[p]);
        migrationCounter[p]++;
    }

    int metapopSize() {
        int tot = Auxils.arraySum(popN);
        // for (int i = 0; i < totSites; i++)
        //     if (alive[i])
        //         tot ++;
        return tot;
    }

    int popSize(int p) {
        int tot = popN[p];
        // for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
        //     if (alive[i])
        //         tot ++;
        return tot;
    }

    double genotypeMean(int t) {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            if (alive[i])
                mean += Auxils.arrayMean(genotype[i], evol.traitGenes[t]);
        mean /= metapopSize();
        return mean;
    }

    double genotypeMean(int p, int t) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                mean += Auxils.arrayMean(genotype[i], evol.traitGenes[t]);
        mean /= popSize(p);
        return mean;
    }

    double genotypeVar(int t) {
        double mean = genotypeMean(t);
        double var = 0;
        for (int i = 0; i < totSites; i++)
            if (alive[i])
                var += Math.pow(mean - Auxils.arrayMean(genotype[i], evol.traitGenes[t]), 2);
        var /= metapopSize();
        return var;
    }

    double genotypeVar(int p, int t) {
        double mean = genotypeMean(p, t);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                var += Math.pow(mean - Auxils.arrayMean(genotype[i], evol.traitGenes[t]), 2);
        var /= popSize(p);
        return var;
    }

    double genotypeMin(int p, int t) {
        double gtp = Auxils.arrayMean(genotype[p * comm.microsites], evol.traitGenes[t]);
        double min = gtp;
        for (int i = p * comm.microsites + 1; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                gtp = Auxils.arrayMean(genotype[i], evol.traitGenes[t]);
                if (gtp < min)
                    min = gtp;
            }
        }
        return min;
    }

    double genotypeMax(int p, int t) {
        double gtp = Auxils.arrayMean(genotype[p * comm.microsites], evol.traitGenes[t]);
        double max = gtp;
        for (int i = p * comm.microsites + 1; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                gtp = Auxils.arrayMean(genotype[i], evol.traitGenes[t]);
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

    double absFitnessMax(int p) {
        double max = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                if (fitness[i] > max)
                    max = fitness[i];
            }
        }
        return max;
    }

    double absContrMean(int p) {
        double mean = 0, relFit = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                relFit = (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
                mean += comm.r*(1 - comm.d)*relFit*(1 - pSex[i]*comm.demogrCost[dcPos]);
            }
        }
        mean /= popSize(p);
        return mean;
    }

    double absContrVar(int p) {
        double mean = absContrMean(p);
        double var = 0, relFit = 0, absContr = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                relFit = (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
                absContr = comm.r*(1 - comm.d)*relFit*(1 - pSex[i]*comm.demogrCost[dcPos]);
                var += Math.pow(mean - absContr, 2);
            }
        }
        var /= popSize(p);
        return var;
    }

    double relContrVar(int p) {
        double mean = absContrMean(p);
        double var = 0, relFit = 0, relContr = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                relFit = (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
                relContr = comm.r*(1 - comm.d)*relFit*(1 - pSex[i]*comm.demogrCost[dcPos])/mean;
                var += Math.pow(1 - relContr, 2);
            }
        }
        var /= popSize(p);
        return var;
    }

    double relRelFitnessVar(int p) {
        double mean = relFitnessMean(p);
        double var = 0, relFit = 0, relRelFit = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                relFit = (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
                relRelFit = relFit/mean;
                var += Math.pow(1 - relRelFit, 2);
            }
        }
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
        int endGens = 0, allInds = popSize(p);
        double nDistinct;
        // int[] allGens = new int[comm.microsites*evol.traitLoci*2];
        int[] allGens = new int[allInds*evol.traitLoci*2];
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            if (alive[i]) {
                System.arraycopy(Auxils.arrayElements(migrationGenotype[i], evol.somGenes), 0, allGens, endGens, 2 * evol.traitLoci);
                endGens += 2 * evol.traitLoci;
            }
        }
        // nDistinct = Auxils.countDistinct(Arrays.copyOf(allGens, endGens));
        nDistinct = Auxils.countDistinct(allGens);
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

    double pDisp() {
        double mean = Auxils.arraySum(pDispSum)/Auxils.arraySum(popN);
        // for (int i = 0; i < totSites; i++)
        //     if (alive[i])
        //         mean += pDisp[i];
        // mean /= metapopSize();
        return mean;
    }

    double pDisp(int p) {
        double mean = pDispSum[p]/popN[p];
        // for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
        //     if (alive[i])
        //         mean += pDisp[i];
        // mean /= popSize(p);
        return mean;
    }

    double pDispVar(int p) {
        double mean = pDisp(p);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (alive[i])
                var += Math.pow(mean - pDisp[i], 2);
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
    // double rho = 1;
    String sexType = "SWITCH";
    double[] pSex = {0.};

    int[] traitDim;

    void init() {
        nbrPatches = gridSize * gridSize;

        traitDim = new int[traits];
        int dim = 0;
        for (int tr = 0; tr < traits; tr++) {
            traitDim[tr] = dim++;
            if (dim == envDims)
                dim = 0;
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
    int dispLoci = 20;
    double dispCorr = 2;
    double dispPow = 4;
    double refDisp = 1e-4;
    int allLoci = traitLoci + sexLoci + dispLoci;
    double mutationRate = 1e-4;
    double mutationRateSex = 1e-5;
    double mutationRateDisp = 1e-4;
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
    int[] dispMother;
    int[] dispFather;
    int[] dispGenes;
    int[] shuffleLocs;
    int[] shuffleTrait;
    int[] shuffleSex;
    int[] shuffleDisp;


    int longPos = 0;

    void init(Comm comm) {
        divF = 2 * Math.pow(Math.sqrt(comm.traits) * omegaE, 2);

        allLoci = traitLoci + sexLoci + dispLoci;

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
        dispMother = new int[dispLoci];
        dispFather = new int[dispLoci];
        dispGenes = new int[2 * dispLoci];
        shuffleLocs = Auxils.enumArray(0, allLoci - 1);
        shuffleTrait = Auxils.enumArray(0, traitLoci*2 - 1);
        shuffleSex = Auxils.enumArray(0, sexLoci*2 - 1);
        shuffleDisp = Auxils.enumArray(0, dispLoci*2 - 1);

        /* somatic genes */
        for (int tr = 0; tr < comm.traits; tr++) {
            for (int l = 0; l < lociPerTrait; l++) {
                longPos = l + (tr * lociPerTrait);
                traitMother[tr][l] = longPos;
                traitFather[tr][l] = traitMother[tr][l] + allLoci;
                somMother[longPos] = longPos;
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

        /* dispersal genes */
        for (int l = 0; l < dispLoci; l++) {
            dispMother[l] = l + traitLoci + sexLoci;
            dispFather[l] = dispMother[l] + allLoci;
        }
        dispGenes = Auxils.arrayConcat(dispMother, dispFather);

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

    public Init(Comm comm, int ps) {
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
        pSex = comm.pSex[ps];
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
                        size = Integer.parseInt(words[1]);
                        comm.pSex = new double[size];
                        for (int i = 0; i < size; i++)
                            comm.pSex[i] = Double.parseDouble(words[2 + i]);
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
                    // case "RHO":
                    //     comm.rho = Double.parseDouble(words[1]);
                    //     break;

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
                    case "DISPLOCI":
                        evol.dispLoci = Integer.parseInt(words[1]);
                        break;
                    case "MUDISP":
                        evol.mutationRateDisp = Double.parseDouble(words[1]);
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
            input.close();
        }
    }
}


/* Auxiliary functions for array calculations */
class Auxils {
//    static UniformRandomProvider random = RandomSource.create(RandomSource.MT_64);
//    static UniformRandomProvider random = RandomSource.create(RandomSource.JSF_64);
//    static UniformRandomProvider random = RandomSource.create(RandomSource.MSWS);
    static UniformRandomProvider random = RandomSource.XO_RO_SHI_RO_128_PP.create();

    static NormalizedGaussianSampler gaussianSampler = ZigguratNormalizedGaussianSampler.of(random);
    static SharedStateDiscreteSampler binomialSamplerSomatic;
    static SharedStateDiscreteSampler binomialSamplerSex;
    static SharedStateDiscreteSampler binomialSamplerDisp;
    static SharedStateDiscreteSampler binomialSamplerInherit;

    static void init(Comm comm, Evol evol) {
        binomialSamplerSomatic = Binomial.of(random, evol.traitLoci*2, evol.mutationRate);
        binomialSamplerSex = Binomial.of(random, evol.sexLoci*2, evol.mutationRateSex);
        binomialSamplerDisp = Binomial.of(random, evol.dispLoci*2, evol.mutationRateDisp);
        binomialSamplerInherit = Binomial.of(random, evol.allLoci, 0.5);
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

    static void arrayShuffle(int[] array, int k) {
        int index, temp;
        for (int i = 0; i < k; i++) {
            index = random.nextInt(i, array.length);
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

    static int[] arraySample(int n, int[] array, int end) {
        int[] tempArr = Arrays.copyOf(array, end);
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

    static int randIntCumProb(double[] cumProbs, int endPos) {
        int val;
        val = Arrays.binarySearch(cumProbs, 0, endPos, random.nextDouble());
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

    static double[] seqArray(double from, double to, double step) {
        int len = (int) Math.round((to - from + step)/step);
        // System.out.println("  from = " + from + "; to = " + to + "; step = " + step + "; len = " + len);
        double[] newArr = new double[len];
        for (int i = 0; i < newArr.length; i++) {
            newArr[i] = from;
            from += step;
        }
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

    static double arrayMean(int[] array, int[] pos) {
        int sum = 0;
        double mean = 0;
        for (int i : pos) sum += array[i];
        mean = ((double) sum)/pos.length;
        return mean;
    }
    
    static double arrayMean(byte[] array, int[] pos) {
        byte sum = 0;
        double mean = 0;
        for (int i : pos) sum += array[i];
        mean = ((double) sum)/pos.length;
        return mean;
    }

    static double arrayMean(double[] array, int[] pos) {
        double mean = 0;
        for (int i : pos) mean += array[i];
        mean /= pos.length;
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

    static void arrayCumSum(double[] array, int endPos) {
        for (int i = 1; i < endPos; i++)
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

    static void arrayDiv(double[] array, int endPos, double a) {
        for (int i = 0; i < endPos; i++)
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
}


