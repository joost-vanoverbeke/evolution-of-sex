

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.CombinationSampler;
import org.apache.commons.rng.sampling.distribution.MarsagliaTsangWangDiscreteSampler.Binomial;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.SharedStateDiscreteSampler;
import org.apache.commons.rng.sampling.distribution.ZigguratNormalizedGaussianSampler;
import org.apache.commons.rng.simple.RandomSource;

import java.io.*;
import java.util.Arrays;

enum OriginType {
    IMM, F1IMM, RES, IMMSEX, RESSEX, MIXIMMSEX, MIXRESSEX
}


/* class EvolvingMetacommunity
 * loops over cycles (time steps) of reproduction and dispersal
 * writes output to file */
public class EvolSex {

    static Comm comm;
    static Evol evol;
    static Run run;

    static Sites sites;

    public static void main(String[] args) throws IOException {

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
                    for (int es = 0; es < comm.envStep.length; es++)
                            for (int dr = 0; dr < comm.dispRate.length; dr++) {

                                System.out.format("run = %d; dims = %d; traits = %d; demCorr = %.2f; disp = %.4f; step = %.4f%n",
                                        (r + 1), comm.envDims, comm.traits, comm.demogrCost[dc], comm.dispRate[dr], comm.envStep[es]);

                                comm.init();
                                evol.init(comm);
                                Auxils.init(comm, evol);
                                Init init = new Init(comm);

                                sites = new Sites(comm, evol, init, dc, es, dr);

                                System.out.format("  time = %d; metacommunity N = %d; absFit = %f; relFit = %f; origin = %f; pSex = %f%n",
                                        0, sites.metaPopSize(), sites.absFitnessMean(), sites.relFitnessMean(), sites.originRatioMean(), sites.pSex());
                                logResults(0, streamOut, r, dc, es, dr);

                                for (int t = 0; t < run.timeSteps; t++) {
                                    sites.changeEnvironment();
                                    sites.findMaxFitness();
                                    sites.contributionAdults();
                                    sites.reproduction(t);

                                    if (t == 0 || ((t + 1) % run.printSteps) == 0) {
                                        System.out.format("  time = %d; metacommunity N = %d; absFit = %f; relFit = %f; origin = %f; pSex = %f%n",
                                                (t + 1), sites.metaPopSize(), sites.absFitnessMean(), sites.relFitnessMean(), sites.originRatioMean(), sites.pSex());
                                    }
                                    if (t == 0 || ((t + 1) % run.saveSteps) == 0) {
                                        sites.findMaxFitness();
                                        logResults(t+1, streamOut, r, dc, es, dr);
                                    }
                                }
                            }

            long endTime = System.currentTimeMillis();
            System.out.println("EvolMetac took " + (endTime - startTime) +
                    " milliseconds.");
        }
    }

    static void logTitles(PrintWriter out) {
        out.print("gridsize;patches;p_e_change;e_step;m;rho;dims;sigma_e;microsites;d;demogr_cost;traits;traitLoci;sigma_z;mu;mu_sex;omega_e;"
//                + "run;time;patch;N;N_res;N_f1imm;N_imm;N_ressex;N_immsex;N_mixressex;N_miximmsex;"
                + "run;time;patch;N;p_migr10;p_migr20;p_migr50;p_migr100;"
//                + "trait_fitness_mean;trait_fitness_var;fitness_mean;fitness_var;fitness_geom;load_mean;load_var;load_geom;S_mean;S_var;pSex_mean;pSex_var;"
                + "pSex_mean;pSex_var;fitness_mean;fitness_var;abs_fitness_mean;load_mean;load_var;S_mean;S_var;"
//                + "origin_mean;origin_max_fitness;origin_sex_mean;origin_sex_max_fitness;origin_asex_mean;origin_asex_max_fitness;"
                + "migration_generation; migration_mother; migration_father; migration_max_fitness; migration_mother_max_fitness; migration_father_max_fitness; migration_genotype; migration_genotype_range;"
                + "pSex_mothers; fitness_mothers; abs_fitness_mothers; migration_mothers;");
//                + "origin_mean;origin_min;origin_max;origin_max_fitness;"
//                + "abs_fitness_mean;abs_fitness_max;abs_fitness_max_res;abs_fitness_max_f1imm;abs_fitness_max_imm;abs_fitness_max_ressex;abs_fitness_max_immsex;abs_fitness_max_mixressex;abs_fitness_max_miximmsex;"
//                + "rel_fitness_max_res;rel_fitness_max_f1imm;rel_fitness_max_imm;rel_fitness_max_ressex;rel_fitness_max_immsex;rel_fitness_max_mixressex;rel_fitness_max_miximmsex");
        for (int tr = 0; tr < comm.traits; tr++)
//            out.format(";dim_tr%d;e_dim_tr%d;genotype_mean_tr%d;genotype_var_tr%d;genotype_min_tr%d;genotype_max_tr%d;phenotype_mean_tr%d;phenotype_var_tr%d;phenotype_min_tr%d;phenotype_max_tr%d;fitness_mean_tr%d;fitness_var_tr%d;"
//                            + "genotype_meta_var_tr%d;phenotype_meta_var_tr%d",
            out.format(";dim_tr%d;e_dim_tr%d;genotype_mean_tr%d;genotype_var_tr%d;phenotype_mean_tr%d;phenotype_var_tr%d;fitness_mean_tr%d;fitness_var_tr%d;"
                            + "genotype_meta_var_tr%d;phenotype_meta_var_tr%d",
                    tr + 1, tr + 1, tr + 1, tr + 1, tr + 1, tr + 1, tr + 1, tr + 1, tr + 1, tr + 1);
        out.println("");
    }

    static void logResults(int t, PrintWriter out, int r, int dc, int es, int dr) {
        for (int p = 0; p < comm.nbrPatches; p++) {
            out.format("%d;%d;%f;%f;%f;%f;%d;%f;%d;%f;%f;%d;%d;%f;%f;%f;%f",
                    comm.gridSize, comm.nbrPatches, comm.pChange, comm.envStep[es], comm.dispRate[dr], comm.rho, comm.envDims, comm.sigmaE, comm.microsites, comm.d, comm.demogrCost[dc], comm.traits, evol.traitLoci, evol.sigmaZ, evol.mutationRate, evol.mutationRateSex, evol.omegaE);
            out.format(";%d;%d;%d;%d;%f;%f;%f;%f",
                    r + 1, t, p + 1, sites.popSize(), sites.pMigr(p, 10), sites.pMigr(p, 20), sites.pMigr(p, 50), sites.pMigr(p, 100));
            out.format(";"
//                            + ";%d;%d;%d;%d;%d;%d;%d;%d;"
//                            + "%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;"
                            + "%f;%f;%f;%f;%f;%f;%f;%f;%f;"
//                            + "%f;%f;%f;%f;%f;%f;"
                            + "%f;%f;%f;%f;%f;%f;%f;%f;"
                            + "%f;%f;%f;%f",
//                            + "%f;%f;%f;%f;"
//                            + "%.10f;%.10f;%.10f;%.10f;%.10f;%.10f;%.10f;%.10f;%.10f;"
//                            + "%.10f;%.10f;%.10f;%.10f;%.10f;%.10f;%.10f",
//                    sites.popSize(), sites.popSizeType(p, OriginType.RES), sites.popSizeType(p, OriginType.F1IMM), sites.popSizeType(p, OriginType.IMM), sites.popSizeType(p, OriginType.RESSEX), sites.popSizeType(p, OriginType.IMMSEX), sites.popSizeType(p, OriginType.MIXRESSEX), sites.popSizeType(p, OriginType.MIXIMMSEX),
//                    sites.traitFitnessMean(p), sites.traitFitnessVar(p), sites.relFitnessMean(p), sites.relFitnessVar(p), sites.relFitnessGeom(p), sites.relLoadMean(p), sites.relLoadVar(p), sites.relLoadGeom(p), sites.selectionDiff(p), sites.selectionDiffVar(p), sites.pSex(p), sites.pSexVar(p),
                    sites.pSex(p), sites.pSexVar(p), sites.relFitnessMean(p), sites.relFitnessVar(p), sites.absFitnessMean(p), sites.relLoadMean(p), sites.relLoadVar(p), sites.selectionDiff(p), sites.selectionDiffVar(p),
//                    sites.originRatioMean(p), sites.originRatioMaxFitness(p), sites.originSexRatioMean(p), sites.originSexRatioMaxFitness(p), sites.originAsexRatioMean(p),sites.originAsexRatioMaxFitness(p),
                    sites.migrationGenerationMean(p), sites.migrationMotherMean(p), sites.migrationFatherMean(p), sites.migrationMaxFitness(p), sites.migrationMotherMaxFitness(p), sites.migrationFatherMaxFitness(p), sites.migrationGenotypeMean(p), sites.migrationGenotypeRange(p),
                    sites.pSexMothers[p], sites.fitnessRelMothers[p], sites.fitnessMothers[p], sites.migrationMothers[p]);
//                    sites.originRatioMean(p),sites.originRatioMin(p),sites.originRatioMax(p),sites.originRatioMaxFitness(p),
//                    sites.absFitnessMean(p), sites.absFitnessMax(p), sites.absFitnessMaxType(p, OriginType.RES), sites.absFitnessMaxType(p, OriginType.F1IMM), sites.absFitnessMaxType(p, OriginType.IMM), sites.absFitnessMaxType(p, OriginType.RESSEX), sites.absFitnessMaxType(p, OriginType.IMMSEX), sites.absFitnessMaxType(p, OriginType.MIXRESSEX), sites.absFitnessMaxType(p, OriginType.MIXIMMSEX),
//                    sites.relFitnessMaxType(p, OriginType.RES), sites.relFitnessMaxType(p, OriginType.F1IMM), sites.relFitnessMaxType(p, OriginType.IMM), sites.relFitnessMaxType(p, OriginType.RESSEX), sites.relFitnessMaxType(p, OriginType.IMMSEX), sites.relFitnessMaxType(p, OriginType.MIXRESSEX), sites.relFitnessMaxType(p, OriginType.MIXIMMSEX));
//            for (int tr = 0; tr < comm.traits; tr++)
//                out.format(";%d;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f",
//                        sites.comm.traitDim[tr] + 1, sites.environment[p][sites.comm.traitDim[tr]], sites.genotypeMean(p, tr), sites.genotypeVar(p, tr), sites.genotypeMin(p, tr), sites.genotypeMax(p, tr), sites.phenotypeMean(p, tr), sites.phenotypeVar(p, tr), sites.phenotypeMin(p, tr), sites.phenotypeMax(p, tr), sites.traitFitnessMean(p, tr), sites.traitFitnessVar(p, tr),
//                        sites.genotypeVar(tr), sites.phenotypeVar(tr));
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
    int[] origin;
    double[] migrationGeneration;
    double[] migrationMother;
    double[] migrationFather;
    double[] migrationSex;
    double[] originRatio;
    double[] originSexRatio;
    double[] originAsexRatio;
    OriginType[] type;
    int[][] newbornsOrigin;
    double[][] newbornsMigrationGeneration;
    double[][] newbornsMigrationMother;
    double[][] newbornsMigrationFather;
    double[][] newbornsMigrationSex;
    double[][] newbornsOriginRatio;
    double[][] newbornsOriginSexRatio;
    double[][] newbornsOriginAsexRatio;
    OriginType[][] newbornsType;
    double[][] traitPhenotype;
    double[][] traitFitness;
    double[] fitness;
    double[] pSex;

    byte[][] genotype;
    byte[][][] newborns;
    int[][] migrationGenotype;
    int[][][] migrationNewborns;

    double[][] environment;
    double[] maxFitness;
    double[] fitnessMothers;
    double[] fitnessRelMothers;
    double[] migrationMothers;
    double[] pSexMothers;

    boolean[] sexAdults;
    int[] endPosFathers;
    int[][] fathersPos;
    double[][] fathersProb;
    double[][] fathersCumProb;
    double[][] mothersCumProb;



    public Sites(Comm cmm, Evol evl, Init init, int dc, int es, int dr) {
        comm = cmm;
        evol = evl;
        dcPos = dc;
        esPos = es;
        drPos = dr;

        comm.calcDispNeighbours(drPos);

        totSites = comm.nbrPatches * comm.microsites;

        patch = new int[totSites];
        origin = new int[totSites];
        migrationGeneration = new double[totSites];
        migrationMother = new double[totSites];
        migrationFather = new double[totSites];
        migrationSex = new double[totSites];
        originRatio = new double[totSites];
        originSexRatio = new double[totSites];
        originAsexRatio = new double[totSites];
        type = new OriginType[totSites];
        traitPhenotype = new double[totSites][comm.traits];
        traitFitness = new double[totSites][comm.traits];
        fitness = new double[totSites];
        genotype = new byte[totSites][2 * evol.allLoci];

        migrationGenotype = new int[totSites][2 * evol.allLoci];

        pSex = new double[totSites];

        newborns = new byte[comm.nbrPatches][comm.nbrNewborns][2 * evol.allLoci];

        migrationNewborns = new int[comm.nbrPatches][comm.nbrNewborns][2 * evol.allLoci];

        newbornsOrigin = new int[comm.nbrPatches][comm.nbrNewborns];
        newbornsMigrationGeneration = new double[comm.nbrPatches][comm.nbrNewborns];
        newbornsMigrationMother = new double[comm.nbrPatches][comm.nbrNewborns];
        newbornsMigrationFather = new double[comm.nbrPatches][comm.nbrNewborns];
        newbornsMigrationSex = new double[comm.nbrPatches][comm.nbrNewborns];
        newbornsOriginRatio = new double[comm.nbrPatches][comm.nbrNewborns];
        newbornsOriginSexRatio = new double[comm.nbrPatches][comm.nbrNewborns];
        newbornsOriginAsexRatio = new double[comm.nbrPatches][comm.nbrNewborns];
        newbornsType = new OriginType[comm.nbrPatches][comm.nbrNewborns];

        environment = new double[comm.nbrPatches][comm.envDims];
        maxFitness = new double[comm.nbrPatches];
        fitnessMothers = new double[comm.nbrPatches];
        fitnessRelMothers = new double[comm.nbrPatches];
        migrationMothers = new double[comm.nbrPatches];
        pSexMothers = new double[comm.nbrPatches];

        sexAdults = new boolean[totSites];
        endPosFathers = new int[comm.nbrPatches];
        fathersPos = new int[comm.nbrPatches][comm.microsites];
        fathersProb = new double[comm.nbrPatches][comm.microsites];
        fathersCumProb = new double[comm.nbrPatches][];
        mothersCumProb = new double[comm.nbrPatches][totSites];

        double indGtp;
        Arrays.fill(maxFitness, 0);
        Arrays.fill(originRatio, 0.);
        Arrays.fill(migrationGeneration, 1.);
        Arrays.fill(migrationMother, 1.);
        Arrays.fill(migrationFather, 1.);
        Arrays.fill(migrationSex, 0.);
        Arrays.fill(originSexRatio, 0.);
        Arrays.fill(originAsexRatio, 0.);

        for (int p = 0; p < comm.nbrPatches; p++) {
            if (comm.envDims >= 0) System.arraycopy(init.environment[p], 0, environment[p], 0, comm.envDims);
            for (int m = (p * comm.microsites); m < ((p + 1) * comm.microsites); m++)
                patch[m] = origin[m] = p;
            int[] posInds = Auxils.arraySample(init.N[p], Auxils.enumArray(p * comm.microsites, ((p + 1) * comm.microsites) - 1));
            for (int m : posInds) {
                type[m] = OriginType.RES;
                fitness[m] = 1;
                for (int tr = 0; tr < comm.traits; tr++) {
                    traitFitness[m][tr] = 1;
                    indGtp = init.genotype[p][tr];
                    for (int l : evol.traitGenes[tr]) {
                        genotype[m][l] = (byte) Math.round(Auxils.random.nextDouble() * 0.5 * (Auxils.random.nextBoolean() ? -1 : 1) + indGtp);
                        migrationGenotype[m][l] = 1;
                    }
                    for (int l : evol.sexGenes) {
                        genotype[m][l] = (byte) ((init.pSex < 0.5) ? 0 : 1);
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

    void changeEnvironment() {
        boolean globalEnv = comm.envType.equals("GLOBAL");
        boolean globalChange;
        double globalStep = 0;
        double step;

        for (int d = 0; d < comm.envDims; d++) {
            globalChange = globalEnv && (Auxils.random.nextDouble() <= comm.pChange);
            if (globalChange) globalStep = comm.envStep[esPos] * (Auxils.random.nextBoolean() ? -1 : 1);
            for (int p = 0; p < comm.nbrPatches; p++) {
                if (globalEnv ? globalChange : (Auxils.random.nextDouble() <= comm.pChange)) {
                    step = globalEnv ? globalStep : (comm.envStep[esPos] * (Auxils.random.nextBoolean() ? -1 : 1));
                    environment[p][d] = environment[p][d] + step;
                    environment[p][d] = Auxils.adjustToRange(environment[p][d], comm.minEnv, comm.maxEnv);
                    adjustFitness(p, d);
                }
            }
        }
    }

    void adjustFitness(int p, int d) {
        double oldFit;
        for (int m = (p * comm.microsites); m < ((p + 1) * comm.microsites); m++) {
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

    void findMaxFitness() {
        Arrays.fill(maxFitness, 0.);
        for (int i = 0; i < totSites; i++)
            if (maxFitness[patch[i]] < fitness[i])
                maxFitness[patch[i]] = fitness[i];
    }

    void contributionAdults() {
        Arrays.fill(endPosFathers, 0);
        double contr;
        int p;

        for (int i = 0; i < totSites; i++) {
            p = patch[i];
            contr = (maxFitness[p] > 0.) ? fitness[i] / maxFitness[p] : 0.;
            sexAdults[i] = Auxils.random.nextDouble() <= pSex[i];
            if (sexAdults[i]) {
                fathersPos[p][endPosFathers[p]] = i;
                fathersProb[p][endPosFathers[p]++] = contr;
                contr *= comm.demogrCost[dcPos];
            }
            for (int p2 = 0; p2 < comm.nbrPatches; p2++) {
                mothersCumProb[p2][i] = contr * comm.dispNeighbours[p2][p];
                if (i > 0)
                    mothersCumProb[p2][i] += mothersCumProb[p2][i - 1];
            }
        }
        for (p = 0; p < comm.nbrPatches; p++) {
            Auxils.arrayDiv(mothersCumProb[p], mothersCumProb[p][totSites-1]);
            if (endPosFathers[p] > 0) {
                fathersCumProb[p] = Arrays.copyOf(fathersProb[p], endPosFathers[p]);
                Auxils.arrayCumSum(fathersCumProb[p]);
                Auxils.arrayDiv(fathersCumProb[p], fathersCumProb[p][endPosFathers[p] - 1]);
            }
        }
    }

    void reproduction(int t) {
        int[] posOffspring;
        int m, f, patchMother, locals;
        Arrays.fill(fitnessMothers, 0.);
        Arrays.fill(fitnessRelMothers, 0.);
        Arrays.fill(migrationMothers, 0.);
        Arrays.fill(pSexMothers, 0.);

        for (int p = 0; p < comm.nbrPatches; p++) {
            locals = 0;
            //sampling parents with replacement!
            //selfing allowed!
            for (int i = 0; i < comm.nbrNewborns; i++) {
                m = Auxils.randIntCumProb(mothersCumProb[p]);
                patchMother = patch[m];
                newbornsOrigin[p][i] = patchMother;
                if (sexAdults[m]) {
                    f = fathersPos[patchMother][Auxils.randIntCumProb(fathersCumProb[patchMother])];
                    if (patchMother != p) {
                        newbornsMigrationGeneration[p][i] = 1.;
                        newbornsMigrationMother[p][i] = 1.;
                        newbornsMigrationFather[p][i] = 1.;
//                        newbornsMigrationGeneration[p][i] = t;
//                        newbornsMigrationMother[p][i] = t;
//                        newbornsMigrationFather[p][i] = t;
//                        newbornsMigrationSex[p][i] = 0.;
//                        newbornsOriginRatio[p][i] = 1.;
//                        newbornsOriginSexRatio[p][i] = 1.;
//                        newbornsOriginAsexRatio[p][i] = 1.;
                    }
                    else {
                        locals++;
                        fitnessMothers[p] += fitness[m];
                        fitnessRelMothers[p] += (maxFitness[p] > 0.) ? fitness[m] / maxFitness[p] : 0.;
                        migrationMothers[p] += Auxils.arrayMean(Auxils.arrayElements(migrationGenotype[m], evol.somGenes));
                        pSexMothers[p] += pSex[m];

                        newbornsMigrationGeneration[p][i] = ((migrationGeneration[m] + migrationGeneration[f]) / 2.) + 1.;
//                        newbornsMigrationGeneration[p][i] = ((migrationGeneration[m] + migrationGeneration[f]) / 2.);
//                        newbornsMigrationSex[p][i] = ((migrationSex[m] + migrationSex[f]) / 2.) + 1.;
//                        newbornsOriginRatio[p][i] = ((originRatio[m] + originRatio[f]) / 2.) * (1 - comm.d);
//                        newbornsOriginSexRatio[p][i] = ((originSexRatio[m] + originSexRatio[f]) / 2.) * (1 - comm.d);
//                        newbornsOriginAsexRatio[p][i] = ((originSexRatio[m] + originSexRatio[f]) / 2. * 0.5) * (1 - comm.d);

//                        newbornsMigrationGeneration[p][i] = Math.min(migrationGeneration[m], migrationGeneration[f]) + 1.;
                        newbornsMigrationMother[p][i] = migrationGeneration[m] + 1.;
                        newbornsMigrationFather[p][i] = migrationGeneration[f] + 1.;
//                        newbornsMigrationMother[p][i] = migrationGeneration[m];
//                        newbornsMigrationFather[p][i] = migrationGeneration[f];
//                        newbornsMigrationSex[p][i] = Math.min(migrationSex[m], migrationSex[f]) + 1.;
//                        newbornsOriginRatio[p][i] = Math.max(originRatio[m], originRatio[f]) * (1. - comm.d);
//                        newbornsOriginSexRatio[p][i] = Math.max(originSexRatio[m], originSexRatio[f]) * (1. - comm.d);
//                        newbornsOriginAsexRatio[p][i] = Math.max(originSexRatio[m], originSexRatio[f]) * (1. - comm.d);

//                        newbornsMigrationGeneration[p][i] = migrationGeneration[m] + 1.;
//                        newbornsMigrationSex[p][i] = migrationSex[m] + 1.;
//                        newbornsOriginRatio[p][i] = originRatio[m] * (1 - comm.d);
//                        newbornsOriginSexRatio[p][i] = originSexRatio[m] * 0.5 * (1 - comm.d);
//                        newbornsOriginAsexRatio[p][i] = originAsexRatio[m] * (1 - comm.d);
                    }
                    if (origin[m] == origin[f]) {
                        if (origin[m] == p)
                            newbornsType[p][i] = OriginType.RESSEX;
                        else
                            newbornsType[p][i] = OriginType.IMMSEX;
                    } else if (origin[m] == p || origin[f] == p)
                        newbornsType[p][i] = OriginType.MIXRESSEX;
                    else
                        newbornsType[p][i] = OriginType.MIXIMMSEX;
                    inherit(p, i, m, f);
                } else {
                    if (patchMother != p) {
                        newbornsMigrationGeneration[p][i] = 1.;
                        newbornsMigrationMother[p][i] = 1.;
                        newbornsMigrationFather[p][i] = 1.;
//                        newbornsMigrationGeneration[p][i] = t;
//                        newbornsMigrationMother[p][i] = t;
//                        newbornsMigrationFather[p][i] = t;
//                        newbornsMigrationSex[p][i] = 0.;
//                        newbornsOriginRatio[p][i] = 1.;
//                        newbornsOriginSexRatio[p][i] = 1.;
//                        newbornsOriginAsexRatio[p][i] = 1.;
                    }
                    else {
                        locals++;
                        fitnessMothers[p] += fitness[m];
                        fitnessRelMothers[p] += (maxFitness[p] > 0.) ? fitness[m] / maxFitness[p] : 0.;
                        migrationMothers[p] += Auxils.arrayMean(Auxils.arrayElements(migrationGenotype[m], evol.somGenes));
                        pSexMothers[p] += pSex[m];

                        newbornsMigrationGeneration[p][i] = migrationGeneration[m] + 1.;
                        newbornsMigrationMother[p][i] = migrationGeneration[m] + 1.;
                        newbornsMigrationFather[p][i] = migrationGeneration[m] + 1.;
//                        newbornsMigrationGeneration[p][i] = migrationGeneration[m];
//                        newbornsMigrationMother[p][i] = migrationGeneration[m];
//                        newbornsMigrationFather[p][i] = migrationGeneration[m];
//                        newbornsOriginRatio[p][i] = originRatio[m] * (1 - comm.d);
//                        newbornsOriginSexRatio[p][i] = originSexRatio[m] * 0.5 * (1 - comm.d);
//                        newbornsOriginAsexRatio[p][i] = originAsexRatio[m] * (1 - comm.d);
                    }
                    if (origin[m] == p)
                        newbornsType[p][i] = OriginType.RES;
                    else
                        newbornsType[p][i] = OriginType.IMM;
                    inherit(p, i, m);
                }
                if (patchMother != p)
                    newbornsType[p][i] = OriginType.F1IMM;
                mutate(p, i);
            }

            fitnessMothers[p] /= locals;
            fitnessRelMothers[p] /= locals;
            migrationMothers[p] /= locals;
            pSexMothers[p] /= locals;
        }
        for (int p = 0; p < comm.nbrPatches; p++) {
            posOffspring = Auxils.combinationSamplerPositionOffspring.sample();
            Auxils.arrayAdd(posOffspring, p * comm.microsites);
            settle(p, posOffspring);
        }
    }

    /* install newborns and inherit traits from the parent(s)
     * including mutation */
    void settle(int p, int[] posOffspring) {
        int pos;
        for (int i = 0; i < comm.nbrNewborns; i++) {
            pos = posOffspring[i];
            origin[pos] = newbornsOrigin[p][i];
            migrationGeneration[pos] = newbornsMigrationGeneration[p][i];
            migrationMother[pos] = newbornsMigrationMother[p][i];
            migrationFather[pos] = newbornsMigrationFather[p][i];
            migrationSex[pos] = newbornsMigrationSex[p][i];
            originRatio[pos] = newbornsOriginRatio[p][i];
            originSexRatio[pos] = newbornsOriginSexRatio[p][i];
            originAsexRatio[pos] = newbornsOriginAsexRatio[p][i];
            type[pos] = newbornsType[p][i];
            System.arraycopy(newborns[p][i], 0, genotype[pos], 0, 2 * evol.allLoci);

            if (origin[pos] != p) {
                for (int l = 0; l < (2*evol.allLoci); l++) {
                    migrationGenotype[pos][l] = 1;
                }
            } else
                System.arraycopy(migrationNewborns[p][i], 0, migrationGenotype[pos], 0, 2 * evol.allLoci);


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
    }

    /* inheritance for asexual reproduction (one parent) */
    void inherit(int p, int posOffspring, int posParent) {
        System.arraycopy(genotype[posParent], 0, newborns[p][posOffspring], 0, 2 * evol.allLoci);

        System.arraycopy(migrationGenotype[posParent], 0, migrationNewborns[p][posOffspring], 0, 2 * evol.allLoci);
        for (int l = 0; l < (2*evol.allLoci); l++) {
            migrationNewborns[p][posOffspring][l] += 1;
        }
    }

    /* inheritance for sexual reproduction (two parent) */
    void inherit(int p, int posOffspring, int posMother, int posFather) {
        for (int l = 0; l < evol.allLoci; l++) {
//            newborns[p][posOffspring][evol.allMother[l]] = genotype[posMother][Auxils.random.nextBoolean() ? evol.allMother[l] : evol.allFather[l]];
//            newborns[p][posOffspring][evol.allFather[l]] = genotype[posFather][Auxils.random.nextBoolean() ? evol.allMother[l] : evol.allFather[l]];

            if(Auxils.random.nextBoolean()) {
                newborns[p][posOffspring][evol.allMother[l]] = genotype[posMother][evol.allMother[l]];
                migrationNewborns[p][posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allMother[l]] + 1;
            } else {
                newborns[p][posOffspring][evol.allMother[l]] = genotype[posMother][evol.allFather[l]];
                migrationNewborns[p][posOffspring][evol.allMother[l]] = migrationGenotype[posMother][evol.allFather[l]] + 1;
            }
            if(Auxils.random.nextBoolean()) {
                newborns[p][posOffspring][evol.allFather[l]] = genotype[posFather][evol.allMother[l]];
                migrationNewborns[p][posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allMother[l]] + 1;
            } else {
                newborns[p][posOffspring][evol.allFather[l]] = genotype[posFather][evol.allFather[l]];
                migrationNewborns[p][posOffspring][evol.allFather[l]] = migrationGenotype[posFather][evol.allFather[l]] + 1;
            }
        }
    }

    void mutate(int p, int posOffspring) {
        int k;
        int[] somMutLocs;

        k = Auxils.binomialSamplerSomatic.sample();
        if (k > 0) {
            CombinationSampler combinationSampler = new CombinationSampler(Auxils.random,evol.traitLoci*2, k);
            somMutLocs = Auxils.arrayElements(evol.somGenes, combinationSampler.sample());
            for (int l : somMutLocs) {
                newborns[p][posOffspring][l] += (Auxils.random.nextBoolean() ? -1 : 1);
            }
        }

        // sex - asex switch
    	if (Auxils.random.nextDouble() <= evol.mutationRateSex)
    	{
            for (int l : evol.sexGenes)
    		{
                newborns[p][posOffspring][l] = (byte) ((newborns[p][posOffspring][l] == 0) ? 1 : 0);
    		}
    	}
    }

    int metaPopSize() {
        return totSites;
    }

    int popSize() {
        return comm.microsites;
    }

    int popSizeType(int p, OriginType t) {
        int N = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (type[i] == t)
                N++;
        return N;
    }

    double genotypeMean(int t) {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            mean += Auxils.arrayMean(Auxils.arrayElements(genotype[i], evol.traitGenes[t]));
        mean /= metaPopSize();
        return mean;
    }

    double genotypeMean(int p, int t) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += Auxils.arrayMean(Auxils.arrayElements(genotype[i], evol.traitGenes[t]));
        mean /= popSize();
        return mean;
    }

    double genotypeVar(int t) {
        double mean = genotypeMean(t);
        double var = 0;
        for (int i = 0; i < totSites; i++)
            var += Math.pow(mean - Auxils.arrayMean(Auxils.arrayElements(genotype[i], evol.traitGenes[t])), 2);
        var /= metaPopSize();
        return var;
    }

    double genotypeVar(int p, int t) {
        double mean = genotypeMean(p, t);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            var += Math.pow(mean - Auxils.arrayMean(Auxils.arrayElements(genotype[i], evol.traitGenes[t])), 2);
        var /= popSize();
        return var;
    }

    double genotypeMin(int p, int t) {
        double gtp = Auxils.arrayMean(Auxils.arrayElements(genotype[p * comm.microsites], evol.traitGenes[t]));
        double min = gtp;
        for (int i = p * comm.microsites + 1; i < (p + 1) * comm.microsites; i++) {
            gtp = Auxils.arrayMean(Auxils.arrayElements(genotype[i], evol.traitGenes[t]));
            if (gtp < min)
                min = gtp;
        }
        return min;
    }

    double genotypeMax(int p, int t) {
        double gtp = Auxils.arrayMean(Auxils.arrayElements(genotype[p * comm.microsites], evol.traitGenes[t]));
        double max = gtp;
        for (int i = p * comm.microsites + 1; i < (p + 1) * comm.microsites; i++) {
            gtp = Auxils.arrayMean(Auxils.arrayElements(genotype[i], evol.traitGenes[t]));
            if (gtp > max)
                max = gtp;
        }
        return max;
    }

    double phenotypeMean(int t) {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            mean += traitPhenotype[i][t];
        mean /= metaPopSize();
        return mean;
    }

    double phenotypeMean(int p, int t) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += traitPhenotype[i][t];
        mean /= popSize();
        return mean;
    }

    double phenotypeVar(int t) {
        double mean = phenotypeMean(t);
        double var = 0;
        for (int i = 0; i < totSites; i++)
            var += Math.pow(mean - traitPhenotype[i][t], 2);
        var /= metaPopSize();
        return var;
    }

    double phenotypeVar(int p, int t) {
        double mean = phenotypeMean(p, t);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            var += Math.pow(mean - traitPhenotype[i][t], 2);
        var /= popSize();
        return var;
    }

    double phenotypeMin(int p, int t) {
        double min = traitPhenotype[p * comm.microsites][t];
        for (int i = p * comm.microsites + 1; i < (p + 1) * comm.microsites; i++) {
            if (traitPhenotype[i][t] < min)
                min = traitPhenotype[i][t];
        }
        return min;
    }

    double phenotypeMax(int p, int t) {
        double max = traitPhenotype[p * comm.microsites][t];
        for (int i = p * comm.microsites + 1; i < (p + 1) * comm.microsites; i++) {
            if (traitPhenotype[i][t] > max)
                max = traitPhenotype[i][t];
        }
        return max;
    }

    double traitFitnessMax(int p, int t) {
        double max = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (max < traitFitness[i][t])
                max = traitFitness[i][t];
        return max;
    }

    double traitFitnessMean(int p, int t) {
        double mean = 0;
        double max = traitFitnessMax(p, t);
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += traitFitness[i][t] / max;
        mean /= popSize();
        return mean;
    }

    double traitFitnessVar(int p, int t) {
        double mean = traitFitnessMean(p, t);
        double max = traitFitnessMax(p, t);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            var += Math.pow(mean - traitFitness[i][t] / max, 2);
        var /= popSize();
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
        var /= popSize();
        return var;
    }

    double absFitnessMean() {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            mean += fitness[i];
        mean /= metaPopSize();
        return mean;
    }

    double absFitnessMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += fitness[i];
        mean /= popSize();
        return mean;
    }

    double absFitnessMax(int p) {
        double max = -1;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (fitness[i] > max)
                max = fitness[i];
        return max;
    }

    double absFitnessMaxType(int p, OriginType t) {
        double max = -1;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (type[i] == t && fitness[i] > max)
                max = fitness[i];
        return max;
    }

    double relFitnessMeanType(int p, OriginType t) {
        double mean = 0;
        double N = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (type[i] == t) {
                mean += (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
                N++;
            }
        if (N > 0)
            mean /= N;
        else
            mean = -1;
        return mean;
    }

    double relFitnessMaxType(int p, OriginType t) {
        double max = -1;
        double relfit;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (type[i] == t) {
                relfit = (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
                if (relfit > max)
                    max = relfit;
            }
        return max;
    }

    double originRatioMean() {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            mean += originRatio[i];
        mean /= metaPopSize();
        return mean;
    }

    double originRatioMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += originRatio[i];
        mean /= popSize();
        return mean;
    }

    double migrationGenerationMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += migrationGeneration[i];
        mean /= popSize();
        return mean;
    }

    double migrationMotherMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += migrationMother[i];
        mean /= popSize();
        return mean;
    }

    double migrationFatherMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += migrationFather[i];
        mean /= popSize();
        return mean;
    }

    double migrationMaxFitness(int p) {
        double migr = 0., n = 0.;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (fitness[i] >= maxFitness[p]) {
                migr += migrationGeneration[i];
                n++;
            }
        migr /= n;
        return migr;
    }

    double migrationMotherMaxFitness(int p) {
        double migr = 0., n = 0.;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (fitness[i] == maxFitness[p]) {
                migr += migrationMother[i];
                n++;
            }
        migr /= n;
        return migr;
    }

    double migrationFatherMaxFitness(int p) {
        double migr = 0., n = 0.;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (fitness[i] == maxFitness[p]) {
                migr += migrationFather[i];
                n++;
            }
        migr /= n;
        return migr;
    }

    double migrationGenotypeMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += Auxils.arrayMean(Auxils.arrayElements(migrationGenotype[i], evol.somGenes));
        mean /= popSize();
        return mean;
    }

    double migrationGenotypeRange(int p) {
        int min = 0, max = 0;
        double range = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            max = 0;
            for (int l : evol.somGenes)
                if (migrationGenotype[i][l] > max)
                    max = migrationGenotype[i][l];
            min = max;
            for (int l : evol.somGenes)
                if (migrationGenotype[i][l] < min)
                    min = migrationGenotype[i][l];
            range += max - min;
        }
        range /= popSize();
        return range;
    }



    double migrationMotherFatherMinMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += Math.min(migrationMother[i], migrationFather[i]);
        mean /= popSize();
        return mean;
    }

    double migrationMotherFatherMeanMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += (migrationMother[i] + migrationFather[i])/2.0;
        mean /= popSize();
        return mean;
    }

    double migrationSexMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += migrationSex[i];
        mean /= popSize();
        return mean;
    }

    double originRatioMax(int p) {
        double max = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (originRatio[i] > max)
                max = originRatio[i];
        return max;
    }

    double originRatioMin(int p) {
        double min = 1;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (originRatio[i] < min)
                min = originRatio[i];
        return min;
    }

    double originRatioMaxFitness(int p) {
        double relfit, fitmax =-1, ratio = -1;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            relfit = (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
            if (relfit > fitmax) {
                fitmax = relfit;
                ratio = originRatio[i];
            }
        }
        return ratio;
    }

    double originSexRatioMean() {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            mean += originSexRatio[i];
        mean /= metaPopSize();
        return mean;
    }

    double originSexRatioMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += originSexRatio[i];
        mean /= popSize();
        return mean;
    }

    double originSexRatioMax(int p) {
        double max = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (originSexRatio[i] > max)
                max = originSexRatio[i];
        return max;
    }

    double originSexRatioMin(int p) {
        double min = 1;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (originSexRatio[i] < min)
                min = originSexRatio[i];
        return min;
    }

    double originSexRatioMaxFitness(int p) {
        double relfit, fitmax =-1, ratio = -1;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            relfit = (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
            if (relfit > fitmax) {
                fitmax = relfit;
                ratio = originSexRatio[i];
            }
        }
        return ratio;
    }

    double originAsexRatioMean() {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            mean += originAsexRatio[i];
        mean /= metaPopSize();
        return mean;
    }

    double originAsexRatioMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += originAsexRatio[i];
        mean /= popSize();
        return mean;
    }

    double originAsexRatioMax(int p) {
        double max = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (originAsexRatio[i] > max)
                max = originAsexRatio[i];
        return max;
    }

    double originAsexRatioMin(int p) {
        double min = 1;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (originAsexRatio[i] < min)
                min = originAsexRatio[i];
        return min;
    }

    double originAsexRatioMaxFitness(int p) {
        double relfit, fitmax =-1, ratio = -1;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            relfit = (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
            if (relfit > fitmax) {
                fitmax = relfit;
                ratio = originAsexRatio[i];
            }
        }
        return ratio;
    }

    double originSexAsexRatioMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += originSexRatio[i]/originAsexRatio[i];
        mean /= popSize();
        return mean;
    }

    double relFitnessMean() {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            mean += (maxFitness[patch[i]] == 0) ? 0 : (fitness[i] / maxFitness[patch[i]]);
        mean /= metaPopSize();
        return mean;
    }

    double relFitnessMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
        mean /= popSize();
        return mean;
    }

    double pMigr(int p, int gen) {
        double n = 0.;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (migrationGeneration[i] <= gen) {
                n++;
            }
        return n/popSize();
    }

    double pMigr(int p, int genStart, int genEnd) {
        double n = 0.;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (migrationGeneration[i] > genStart && migrationGeneration[i] <= genEnd) {
                n++;
            }
        return n/popSize();
    }

    double pRes(int p) {
        double n = 0.;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (migrationGeneration[i] > 10) {
                n++;
            }
        return n/popSize();
    }

    double relFitnessMigrMean(int p, int gen) {
        double mean = 0, n = 0.;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (migrationGeneration[i] <= gen) {
                mean += (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
                n++;
            }
        mean /= n;
        return mean;
    }

    double relFitnessMigrMean(int p, int genStart, int genEnd) {
        double mean = 0, n = 0.;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (migrationGeneration[i] > genStart && migrationGeneration[i] <= genEnd) {
                mean += (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
                n++;
            }
        mean /= n;
        return mean;
    }

    double relFitnessResMean(int p, int gen) {
        double mean = 0, n = 0.;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            if (migrationGeneration[i] > gen) {
                mean += (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
                n++;
            }
        mean /= n;
        return mean;
    }

    double relFitnessGeom(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += Math.log((maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]));
        mean /= popSize();
        return Math.exp(mean);
    }

    double relFitnessVar(int p) {
        double mean = relFitnessMean(p);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            var += (maxFitness[p] == 0) ? 0 : Math.pow(mean - fitness[i] / maxFitness[p], 2);
        var /= popSize();
        return var;
    }

    double relFitnessMin(int p) {
        double relFit, min = 1;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            relFit = (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
            if (relFit < min)
                min = relFit;
        }
        return min;
    }

    double relFitnessMax(int p) {
        double relFit, max = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
            relFit = (maxFitness[p] == 0) ? 0 : (fitness[i] / maxFitness[p]);
            if (relFit > max)
                max = relFit;
        }
        return max;
    }

    double relLoadMean(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += (maxFitness[p] == 0) ? 0 : (1 - fitness[i] / maxFitness[p]);
        mean /= popSize();
        return mean;
    }

    double relLoadGeom(int p) {
        double mean = 0;
        double min = 1;
        if (maxFitness[p] == 0)
            mean = 1;
        else {
            for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
                if ((fitness[i] < maxFitness[p]) && ((1 - fitness[i] / maxFitness[p]) < min))
                    min = 1 - fitness[i] / maxFitness[p];
            min = Math.exp(Math.floor(Math.log(min)));
            for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
                mean += Math.log(1 - fitness[i] / maxFitness[p] + min);
            mean /= popSize();
            mean = Math.exp(mean) - min;
        }
        return mean;
    }

    double relLoadVar(int p) {
        double mean = relLoadMean(p);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            var += (maxFitness[p] == 0) ? 0 : Math.pow(mean - (1 - fitness[i] / maxFitness[p]), 2);
        var /= popSize();
        return var;
    }

    double selectionDiff(int p, int t) {
        double mean = 0;
        double sum = 0;
        double fitRel;
        double fitMean = relFitnessMean(p);
        double phenotpSd = Math.sqrt(phenotypeVar(p, t));
        double SDiff;

        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) {
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
        var /= popSize();
        return var;
    }

    double pSex() {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
            mean += pSex[i];
        mean /= metaPopSize();
        return mean;
    }

    double pSex(int p) {
        double mean = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            mean += pSex[i];
        mean /= popSize();
        return mean;
    }

    double pSexVar(int p) {
        double mean = pSex(p);
        double var = 0;
        for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++)
            var += Math.pow(mean - pSex[i], 2);
        var /= popSize();
        return var;
    }
}


/* Ecological parameters/variables */
class Comm {
    String envType = "GLOBAL";
    int envDims = 1;
    int traits = 2;
    double minEnv = 0.2;
    double maxEnv = 0.8;
    double sigmaE = 0.0;
    int microsites = 600;
    double d = 0.1;
    int nbrNewborns = (int) Math.round(microsites*d);
    double[] demogrCost = {0.5};

    int gridSize = 2;
    int nbrPatches = gridSize * gridSize;
    double pChange = 0.1;
    double[] envStep = {0.01};
    double[] dispRate = {0.01};
    double rho = 1;
    double pSex = 0;

    int[] traitDim;

    double[][] neighbours = new double[nbrPatches][nbrPatches];
    double[][] dispNeighbours = new double[nbrPatches][nbrPatches];

    void init() {
        nbrNewborns = (int) (microsites*d);
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

        if (comm.envType.equals("GLOBAL")) {
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
        pSex = comm.pSex >= 0 ? comm.pSex : Auxils.random.nextDouble();
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
                    case "PCHANGE":
                        comm.pChange = Double.parseDouble(words[1]);
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
                    case "SEXLOCI":
                        evol.sexLoci = Integer.parseInt(words[1]);
                        break;
                    case "MU":
                        evol.mutationRate = Double.parseDouble(words[1]);
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
    static CombinationSampler combinationSamplerPositionOffspring;

    static void init(Comm comm, Evol evol) {
        binomialSamplerSomatic = Binomial.of(random, evol.traitLoci*2, evol.mutationRate);
        combinationSamplerPositionOffspring = new CombinationSampler(random, comm.microsites, comm.nbrNewborns);
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

    static int countDistinct(byte[] arr, int n) {
        // First sort the array so that all
        // occurrences become consecutive
        Arrays.sort(arr);

        // Traverse the sorted array
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


