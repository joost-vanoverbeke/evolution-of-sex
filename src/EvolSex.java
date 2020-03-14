

//import java.util.Arrays;

import java.util.Random;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.IOException;


/* class EvolvingMetacommunity
 * loops over cycles (time steps) of patch extinction, mortality, reproduction and dispersal
 * writes output to file */
public class EvolSex {

	public static void main(String[] args) throws IOException {
		
		Comm comm = new Comm();
		Evol evol = new Evol();
		Run run = new Run();
		if (args.length > 0)
			Reader.readInput(args[0], comm, evol, run);

		try (PrintWriter streamOut = new PrintWriter(new FileWriter(run.fileName))) {

			long startTime = System.currentTimeMillis();

			streamOut.print("gridsize;patches;p_ext;p_e_change;e_step;e_dim_dev;e_patch_dev;m;rho;dims;sigma_e;microsites;d;demogr_cost;traits;traitLoci;sigma_z;mu;omega_e;" 
					+ "run;time;patch;N;trait_fitness_mean;trait_fitness_var;fitness_mean;fitness_var;load_mean;load_var;S_mean;S_var;pSex_mean;pSex_var");

			for (int tr = 0; tr < comm.traits; tr++)
        		streamOut.format(";dim_tr%d;e_dim_tr%d;genotype_mean_tr%d;genotype_var_tr%d;phenotype_mean_tr%d;phenotype_var_tr%d;alleles_locus_tr%d;alleles_trait_tr%d;"
        				+ "genotype_meta_var_tr%d;phenotype_meta_var_tr%d;alleles_locus_meta_tr%d;alleles_trait_meta_tr%d",
        				tr+1,tr+1,tr+1,tr+1,tr+1,tr+1,tr+1,tr+1,tr+1,tr+1,tr+1,tr+1);
        	streamOut.println("");
        	
			for (int r = 0; r < run.runs; r++)
			for (int dc = 0; dc < comm.demogrCost.length; dc++) 
			for (int es = 0; es < comm.envStep.length; es++) 
			for (int pe = 0; pe < comm.pExt.length; pe++) 
			for (int dd = 0; dd < comm.dimDev.length; dd++) 
			for (int pd = 0; pd < comm.patchDev.length; pd++) 
			for (int dr = 0; dr < comm.dispRate.length; dr++) 
			{

				System.out.format("run = %d; dims = %d; traits = %d; demCorr = %.2f; step = %.4f; pExt = %.4f; dDev = %.4f; pDev = %.4f; disp = %.4f%n", 
						(r+1), comm.envDims, comm.traits, comm.demogrCost[dc], comm.envStep[es], comm.pExt[pe], comm.dimDev[dd], comm.patchDev[pd], comm.dispRate[dr]);

				comm.init();
				evol.init(comm);
				Init init = new Init(comm, evol, dd, pd, es);

				Sites sites = new Sites(comm, evol, init, dc, es, pe, dd, pd, dr);

				System.out.format("  time = %d; metacommunity N = %d; absFit = %f; fit = %f; pSex = %f%n", 0, sites.metaPopSize(), sites.absFitness(), sites.fitness(), sites.pSex());
				for (int p = 0; p < comm.nbrPatches; p++)
				{
					streamOut.format("%d;%d;%f;%f;%f;%f;%f;%f;%f;%d;%f;%d;%f;%f;%d;%d;%f;%f;%f",
							comm.gridSize, comm.nbrPatches, comm.pExt[pe], comm.pChange, comm.envStep[es], comm.dimDev[dd], comm.patchDev[pd], comm.dispRate[dr], comm.rho, comm.envDims, comm.sigmaE, comm.microsites, comm.d, comm.demogrCost[dc], comm.traits, evol.traitLoci, evol.sigmaZ, evol.mutationRate, evol.omegaE);
					streamOut.format(";%d;%d;%d",
							r + 1, 0, p + 1);
					streamOut.format(";%d;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f",
							sites.popSize(p), sites.traitFitness(p), sites.traitFitnessVar(p), sites.fitness(p), sites.fitnessVar(p), sites.load(p), sites.loadVar(p), sites.selectionDiff(p), sites.selectionDiffVar(p), sites.pSex(p), sites.pSexVar(p));
					for (int tr = 0; tr < comm.traits; tr++)
						streamOut.format(";%d;%f;%f;%f;%f;%f;%d;%d;%f;%f;%d;%d",
								sites.comm.traitDim[tr] + 1, sites.patchDimEnv[p][sites.comm.traitDim[tr]], sites.genotype(p, tr), sites.genotypeVar(p, tr), sites.phenotype(p, tr), sites.phenotypeVar(p, tr), sites.distinctAllelesLocus(p, tr), sites.distinctAllelesTrait(p, tr),
								sites.genotypeVar(tr), sites.phenotypeVar(tr), sites.distinctAllelesLocus(tr), sites.distinctAllelesTrait(tr));
					streamOut.println("");
				}

				for (int t = 0; t < run.timeSteps; t++) {
					sites.patchExtinction();
					sites.changeEnvironment();
					sites.mortality();
			    	sites.contributionAdults();
			        sites.reproduction();

					if (t == 0 || ((t+1) % run.printSteps) == 0) 
					{
						System.out.format("  time = %d; metacommunity N = %d; absFit = %f; fit = %f; pSex = %f%n", (t+1), sites.metaPopSize(), sites.absFitness(), sites.fitness(), sites.pSex());
					}
					if (t == 0 || ((t+1) % run.saveSteps) == 0) 
					{
						for (int p = 0; p < comm.nbrPatches; p++) 
						{
							streamOut.format("%d;%d;%f;%f;%f;%f;%f;%f;%f;%d;%f;%d;%f;%f;%d;%d;%f;%f;%f",
									comm.gridSize,comm.nbrPatches,comm.pExt[pe],comm.pChange,comm.envStep[es],comm.dimDev[dd],comm.patchDev[pd],comm.dispRate[dr],comm.rho,comm.envDims,comm.sigmaE,comm.microsites,comm.d,comm.demogrCost[dc],comm.traits,evol.traitLoci,evol.sigmaZ,evol.mutationRate,evol.omegaE);
							streamOut.format(";%d;%d;%d",
									r+1,t+1,p+1);
							streamOut.format(";%d;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f",
									sites.popSize(p),sites.traitFitness(p),sites.traitFitnessVar(p),sites.fitness(p),sites.fitnessVar(p),sites.load(p),sites.loadVar(p),sites.selectionDiff(p),sites.selectionDiffVar(p),sites.pSex(p),sites.pSexVar(p));
							for (int tr = 0; tr < comm.traits; tr++)
								streamOut.format(";%d;%f;%f;%f;%f;%f;%d;%d;%f;%f;%d;%d",
										sites.comm.traitDim[tr]+1,sites.patchDimEnv[p][sites.comm.traitDim[tr]],sites.genotype(p,tr),sites.genotypeVar(p,tr),sites.phenotype(p,tr),sites.phenotypeVar(p,tr),sites.distinctAllelesLocus(p,tr),sites.distinctAllelesTrait(p,tr),
										sites.genotypeVar(tr),sites.phenotypeVar(tr),sites.distinctAllelesLocus(tr),sites.distinctAllelesTrait(tr));
							streamOut.println("");
						}
					}
				}
			}

			long endTime = System.currentTimeMillis();
			System.out.println("EvolMetac took " + (endTime - startTime) +
					" milliseconds.");

		}
	}
}


/* class Sites
 * keeps track of individuals and their attributes in microsites (microhabitats within patches)
 * implements patch extinction, mortality, reproduction (with inheritance and mutation) and dispersal */
class Sites 
{
    Comm comm;
    Evol evol;
    int totSites;

    int dcPos;
    int esPos;
    int pePos;
    int drPos;

    int[] patch;
    double[][] environment;
    boolean[] alive;
    double[][] traitFenotype;
    double[][] traitFitness;
    double[] fitness;
    double[][] genotype;
    double[] pSex;

    int[] posAdults;
    boolean[] sexAdults;
    int[] nbrAdults;
    int[] cumsumAdults;
    int endPosAdults;
    int[][] posEmpty;
    int[] nbrEmpty;
    double[] maxFitness;

    double[][] adultsProb;

    double[][] patchDimEnv;

    int nbrSettled;
    int[] posOffspring, posMothers;
    boolean[] sexMothers;

    
    public Sites (Comm cmm, Evol evl, Init init, int dc, int es, int pe, int dd, int pd, int dr) 
    {
        comm = cmm;
        evol = evl;
        dcPos = dc;
        esPos = es;
        pePos = pe;
        drPos = dr;

        comm.calcDispNeighbours(drPos);

        totSites = comm.nbrPatches * comm.microsites;
        
        patch = new int [totSites];
        environment = new double [totSites][comm.envDims];
        alive = new boolean [totSites];
        traitFenotype = new double [totSites][comm.traits];
        traitFitness = new double [totSites][comm.traits];
        fitness = new double [totSites];
        genotype = new double [totSites][2*evol.allLoci];
        pSex = new double [totSites];

        
        posAdults = new int[totSites];
        sexAdults = new boolean[totSites];
        nbrAdults = new int[comm.nbrPatches];
        cumsumAdults = new int[comm.nbrPatches];
        posEmpty = new int[comm.nbrPatches][comm.microsites];
        nbrEmpty = new int[comm.nbrPatches];
        maxFitness = new double[comm.nbrPatches];

        adultsProb = new double[comm.nbrPatches][totSites];
        patchDimEnv = new double[comm.nbrPatches][comm.envDims];
        
        double indGtp;
		java.util.Arrays.fill(alive,false);
		java.util.Arrays.fill(maxFitness,0);

        for (int p = 0; p < comm.nbrPatches; p++) {
			if (comm.envDims >= 0) System.arraycopy(init.environment[p], 0, patchDimEnv[p], 0, comm.envDims);
            for (int m = (p * comm.microsites); m < ((p + 1) * comm.microsites); m++) {
                patch[m] = p;
            	for (int d = 0; d < comm.envDims; d++)
            		environment[m][d] = patchDimEnv[p][d] + (Auxils.random.nextGaussian() * comm.sigmaE);
            }
            int[] posInds = Auxils.arraySample(init.N[p], Auxils.enumArray(p * comm.microsites, ((p + 1) * comm.microsites) - 1));
	    	for (int m : posInds) {
            	alive[m] = true;
            	fitness[m] = 1;
            	for (int tr = 0; tr < comm.traits; tr++) {
            		traitFitness[m][tr] = 1;
            		indGtp = init.genotype[p][tr];
            		for (int l : evol.traitGenes[tr]) {
            			genotype[m][l] = (int) Math.round(Auxils.random.nextDouble()*0.5*(Auxils.random.nextBoolean() ? -1 : 1) + indGtp);
            		}
                	for (int l : evol.sexGenes) {
                		genotype[m][l] = (int) Math.round(Auxils.random.nextDouble()*0.5*(Auxils.random.nextBoolean() ? -1 : 1) + init.pSex);
// sex - asex switch
//                		genotype[m][l] =  init.pSex < 0.5 ? 0 : 1;
                	}
            		traitFenotype[m][tr] = Auxils.arrayMean(Auxils.arrayElements(genotype[m],evol.traitGenes[tr])) + (Auxils.random.nextGaussian() * evol.sigmaZ);
            		traitFitness[m][tr] = Math.exp(-(Math.pow(traitFenotype[m][tr]-environment[m][comm.traitDim[tr]],2))/evol.divF);
            		fitness[m] *= traitFitness[m][tr];
            	}
            	if (maxFitness[p] < fitness[m])
            		maxFitness[p] = fitness[m];

        		pSex[m] = Math.min(1, Math.max(0, Auxils.arrayMean(Auxils.arrayElements(genotype[m],evol.sexGenes))));
	    	}
        }
        
//debug
//      for (int m = 0; m < totSites; m++) {
//    	  for (int t = 0; t < comm.traits; t++)
//    		  System.out.println("traitGenes " + m + " " + t + " = " + java.util.Arrays.toString(Auxils.arrayElements(genotype[m],evol.traitGenes[t])));
//    	  System.out.println("sexGenes " + m + " = " + java.util.Arrays.toString(Auxils.arrayElements(genotype[m],evol.sexGenes)));
//      }
//
    }
    
    void patchExtinction()
    {
        for (int p = 0; p < comm.nbrPatches; p++)
            if (Auxils.random.nextDouble() <= comm.pExt[pePos])
                for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
                    alive[i] = false;
    }
    
    void changeEnvironment()
    {
    	double step;
    	if (comm.envSync.equals("YES") && Auxils.random.nextDouble() <= comm.pChange)
    	{
			for (int d = 0; d < comm.envDims; d++)
			{
				if (comm.envType.equals("GLOBAL")) 
				{
					step = comm.envStep[esPos]*(Auxils.random.nextBoolean() ? -1 : 1);
					//    		   			step = comm.envStep[esPos]*Auxils.random.nextGaussian();
					for (int p = 0; p < comm.nbrPatches; p++) 
					{
						patchDimEnv[p][d] = patchDimEnv[p][d] + step;
						patchDimEnv[p][d] = Auxils.adjustToRange(patchDimEnv[p][d], 
								comm.minEnv, 
								comm.maxEnv);
						changeMicrositeEnv(p,d);
					}
				}
				else 
				{
					for (int p = 0; p < comm.nbrPatches; p++)
					{
						step = comm.envStep[esPos]*(Auxils.random.nextBoolean() ? -1 : 1);
						//    		   			step = comm.envStep[esPos]*Auxils.random.nextGaussian();
						patchDimEnv[p][d] = patchDimEnv[p][d] + step;
						patchDimEnv[p][d] = Auxils.adjustToRange(patchDimEnv[p][d], 
								comm.minEnv, 
								comm.maxEnv);
						changeMicrositeEnv(p,d);
					}
				}
			}
    	}
    	else
    	{
    		for (int d = 0; d < comm.envDims; d++)
    		{
    			if (comm.envType.equals("GLOBAL")) 
    			{
    				if (Auxils.random.nextDouble() <= comm.pChange) 
    				{
    					step = comm.envStep[esPos]*(Auxils.random.nextBoolean() ? -1 : 1);
//    		   			step = comm.envStep[esPos]*Auxils.random.nextGaussian();
    					for (int p = 0; p < comm.nbrPatches; p++) 
    					{
    						patchDimEnv[p][d] = patchDimEnv[p][d] + step;
    						patchDimEnv[p][d] = Auxils.adjustToRange(patchDimEnv[p][d], 
    								comm.minEnv, 
    								comm.maxEnv);
    						changeMicrositeEnv(p,d);
    					}
    				}
    			}
    			else 
    			{
    				for (int p = 0; p < comm.nbrPatches; p++)
    				{
    					if (Auxils.random.nextDouble() <= comm.pChange) 
    					{
    						step = comm.envStep[esPos]*(Auxils.random.nextBoolean() ? -1 : 1);
//      		   			step = comm.envStep[esPos]*Auxils.random.nextGaussian();
    						patchDimEnv[p][d] = patchDimEnv[p][d] + step;
    						patchDimEnv[p][d] = Auxils.adjustToRange(patchDimEnv[p][d], 
    								comm.minEnv, 
    								comm.maxEnv);
    						changeMicrositeEnv(p,d);
    					}
    				}
    			}
    		}
    	}
    }

    void changeMicrositeEnv(int p, int d)
    {
		for (int m = (p * comm.microsites); m < ((p + 1) * comm.microsites); m++) 
		{
			environment[m][d] = patchDimEnv[p][d] + (Auxils.random.nextGaussian() * comm.sigmaE);
   			if (alive[m]) 
			{
   				fitness[m] = 1;
   				for (int tr = 0; tr < comm.traits; tr++)
   				{
   					if (comm.traitDim[tr] == d) {
   						traitFitness[m][tr] = Math.exp(-(Math.pow(traitFenotype[m][tr]-environment[m][d],2))/evol.divF);
   					}
   					fitness[m] *= traitFitness[m][tr];
   				}
			}
		}
    }
    

/* mortality 
 * keeps also track of surviving individuals and empty microsites per patch for efficiency in further calculations of reproduction and dispersal */    
    void mortality()
    {

        double maxFitnessTemp;
        endPosAdults = 0;

      java.util.Arrays.fill(nbrAdults, 0);
      java.util.Arrays.fill(nbrEmpty, 0);

    	for (int p = 0; p < comm.nbrPatches; p++) 
    	{
            maxFitnessTemp = 0;
            
            nbrAdults[p] = 0;
            nbrEmpty[p] = 0;
            
            for (int i = p * comm.microsites; i < (p + 1) * comm.microsites; i++) 
            {
            	if (alive[i])
            	{
            		alive[i] = Auxils.random.nextDouble() < (1-comm.d);
            	}
                if (alive[i]) 
                {
                	if (maxFitnessTemp < fitness[i])
                		maxFitnessTemp = fitness[i];
                    posAdults[endPosAdults++] = i;
                    nbrAdults[p]++;
                }
                else 
                {
                	posEmpty[p][nbrEmpty[p]++] = i;
                }
            }
            maxFitness[p] = maxFitnessTemp;
    	}
    	
        System.arraycopy(nbrAdults, 0, cumsumAdults, 0, nbrAdults.length);
        Auxils.arrayCumSum(cumsumAdults);

    }

    void contributionAdults()
    {
    	double contr;
    	for (int p = 0; p < comm.nbrPatches; p++) 
    	{
    		for (int i = (p == 0) ? 0 : cumsumAdults[p-1]; i < cumsumAdults[p]; i++)
    		{
    			sexAdults[i] = Auxils.random.nextDouble() <= pSex[posAdults[i]];
    			if(maxFitness[p] > 0.)
    				contr = (fitness[posAdults[i]]/maxFitness[p]);
    			else
//    				contr = 0.;
    				contr = 1.;

    			if (sexAdults[i])
    				contr *= comm.demogrCost[dcPos];
    			for (int p2 = 0; p2 < comm.nbrPatches; p2++)
    			{
    				adultsProb[p2][i] = contr*comm.dispNeighbours[p2][p];
    			}
    		}
    	}
    }

    void reproduction()
    {   

        int[] sampleM;

    	for (int p = 0; p < comm.nbrPatches; p++) 
    	{
    		nbrSettled = nbrEmpty[p];
    		if (nbrSettled > 0) {
    			posOffspring = Auxils.arraySample(nbrSettled, java.util.Arrays.copyOf(posEmpty[p], nbrEmpty[p]));
    		    //sampling mothers with replacement!
				sampleM = Auxils.arraySampleProb(nbrSettled, Auxils.enumArray(0, endPosAdults-1), java.util.Arrays.copyOf(adultsProb[p], endPosAdults), true);
    			posMothers = Auxils.arrayElements(posAdults, sampleM);
                sexMothers = Auxils.arrayElements(sexAdults, sampleM);
                settle(p);
    		}
    	}
    }    

    /* installs newborns in empty microsites and inherits traits from the parent(s)
     * including mutation */
    void settle(int p)
    {
        int startPosSample, rangeSample;
    	int posFather, nbrFathers, patchMother;
    	int[] posFathers;
    	double[] FathersProb;
    	for (int i = 0; i < nbrSettled; i++) {
    		alive[posOffspring[i]] = true;
    		fitness[posOffspring[i]] = 1;
    		if (sexMothers[i]) {
    			patchMother = patch[posMothers[i]];
    			startPosSample = (patchMother == 0) ? 0 : cumsumAdults[patchMother-1];
    			rangeSample = cumsumAdults[patchMother] - startPosSample; 
    			posFathers = new int[rangeSample];
    			FathersProb = new double[rangeSample];
    			nbrFathers = 0;
    			// selfing allowed
    			for (int f = startPosSample; f < (startPosSample + rangeSample); f++)
    				if (sexAdults[f]) {
    					posFathers[nbrFathers] = posAdults[f];
    					//        					FathersProb[nbrFathers++] = adultsProb[f];
    					FathersProb[nbrFathers++] = adultsProb[p][f];
    				}
    			posFather = posFathers[Auxils.randIntProb(nbrFathers, FathersProb)];
    			inherit(posOffspring[i], posMothers[i], posFather);
    			// no selfing allowed
    			//        			for (int f = startPosSample; f < (startPosSample + rangeSample); f++)
    			//        				if (sexAdults[f] && ID[posAdults[f]] == ID[posMothers[i]] && posAdults[f] != posMothers[i]) {
    			//        					posFathers[nbrFathers] = posAdults[f];
    			//        					FathersProb[nbrFathers++] = adultsProb[f];
    			//        				}
    			//        			if (nbrFathers > 0) {
    			//        				posFather = posFathers[Auxils.randIntProb(nbrFathers, FathersProb)];
    			//        				inherit(posOffspring[i], posMothers[i], posFather);
    			//        			}
    			//        			else 
    			//        				inherit(posOffspring[i], posMothers[i]);
    		}
    		else
    			inherit(posOffspring[i], posMothers[i]);

    		mutate(posOffspring[i]);

    		for (int tr = 0; tr < comm.traits; tr++) {
    			traitFenotype[posOffspring[i]][tr] = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring[i]],evol.traitGenes[tr])) + (Auxils.random.nextGaussian() * evol.sigmaZ);
    			traitFitness[posOffspring[i]][tr] = Math.exp(-(Math.pow(traitFenotype[posOffspring[i]][tr]-environment[posOffspring[i]][comm.traitDim[tr]],2))/evol.divF);
    			fitness[posOffspring[i]] *= traitFitness[posOffspring[i]][tr];
    		}
    		if (maxFitness[p] < fitness[posOffspring[i]])
    			maxFitness[p] = fitness[posOffspring[i]];
    		
    		pSex[posOffspring[i]] = Math.min(1, Math.max(0, Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring[i]],evol.sexGenes))));
    	}
    }

    /* inheritance for asexual reproduction (one parent) */    
    void inherit(int posOffspring, int posParent)
    {
    	System.arraycopy(genotype[posParent], 0, genotype[posOffspring], 0, 2*evol.allLoci);
    }

    /* inheritance for sexual reproduction (two parent) */    
    void inherit(int posOffspring, int posMother, int posFather)
    {
    	for (int l = 0; l < evol.allLoci; l++) {
    		if (Auxils.random.nextBoolean())
    			genotype[posOffspring][evol.allMother[l]] = genotype[posMother][evol.allMother[l]];
    		else
    			genotype[posOffspring][evol.allMother[l]] = genotype[posMother][evol.allFather[l]];
    		if (Auxils.random.nextBoolean())
    			genotype[posOffspring][evol.allFather[l]] = genotype[posFather][evol.allMother[l]];
    		else
    			genotype[posOffspring][evol.allFather[l]] = genotype[posFather][evol.allFather[l]];
    	}
}

    void mutate(int posOffspring) {

    	for (int l : evol.somGenes) {
    		if (Auxils.random.nextDouble() <= evol.mutationRate)
    			genotype[posOffspring][l] += (Auxils.random.nextBoolean() ? -1 : 1)*evol.mutationSize;
//    			genotype[posOffspring][l] += Auxils.random.nextGaussian()*0.5;
    	}

    	double pSexTemp = Auxils.arrayMean(Auxils.arrayElements(genotype[posOffspring],evol.sexGenes));
    	for (int l : evol.sexGenes) {
    		if (Auxils.random.nextDouble() <= evol.mutationRate) {
    			if (pSexTemp <= 0.)
    				genotype[posOffspring][l] += 1;
    			else if (pSexTemp >= 1.)
    				genotype[posOffspring][l] -= 1;
    			else
    				genotype[posOffspring][l] += (Auxils.random.nextBoolean() ? -1 : 1);
//    				genotype[posOffspring][l] += Auxils.random.nextGaussian()*0.5;
    		}
    	}

// sex - asex switch    	
//    	if (Auxils.random.nextDouble() <= evol.mutationRate) 
//    	{
//    		for (int l : evol.sexGenes) 
//    		{
//    			genotype[posOffspring][l] = (genotype[posOffspring][l] == 0) ? 1 : 0;
//    		}
//    	}

    }

    int metaPopSize()
    {
        return Auxils.arraySum(alive);
    }

    int popSize(int p)
    {
		return Auxils.arraySum(java.util.Arrays.copyOfRange(alive, p*comm.microsites, p*comm.microsites + comm.microsites));
    }

    int distinctAllelesLocus(int t)
    {
        int n = 0;
        int end = 0;
        double[] alleles = new double [totSites*2*evol.lociPerTrait];
        for (int l = 0; l < evol.lociPerTrait; l++)
        {
        	for (int i = 0; i < totSites; i++)
        		if (alive[i])
        		{
        			alleles[end++] = genotype[i][evol.traitMother[t][l]]; 
        			alleles[end++] = genotype[i][evol.traitFather[t][l]];
        		}
        	n += Auxils.countDistinct(alleles, end);
        }
        return n;
    }

    int distinctAllelesLocus(int p, int t)
    {
        int n = 0;
        int end = 0;
        double[] alleles = new double [comm.microsites*2*evol.lociPerTrait];
        for (int l = 0; l < evol.lociPerTrait; l++)
        {
        	for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        		if (alive[i])
        		{
        			alleles[end++] = genotype[i][evol.traitMother[t][l]]; 
        			alleles[end++] = genotype[i][evol.traitFather[t][l]];
        		}
        	n += Auxils.countDistinct(alleles, end);
        }
        return n;
    }

    int distinctAllelesTrait(int t)
    {
        int n;
        int end = 0;
        double[] alleles = new double [totSites*2*evol.lociPerTrait];
        for (int i = 0; i < totSites; i++)
        	if (alive[i])
        	{
        		System.arraycopy(Auxils.arrayElements(genotype[i],evol.traitGenes[t]), 0, alleles, end, 2*evol.lociPerTrait);
        		end += 2*evol.lociPerTrait;
        	}
        n = Auxils.countDistinct(alleles, end);
        return n;
    }

    int distinctAllelesTrait(int p, int t)
    {
        int n;
        int end = 0;
        double[] alleles = new double [comm.microsites*2*evol.lociPerTrait];
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        	{
        		System.arraycopy(Auxils.arrayElements(genotype[i],evol.traitGenes[t]), 0, alleles, end, 2*evol.lociPerTrait);
        		end += 2*evol.lociPerTrait;
        	}
        n = Auxils.countDistinct(alleles, alleles.length);
        return n;
    }

    double genotype(int t)
    {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
        	if (alive[i])
        		mean += Auxils.arrayMean(Auxils.arrayElements(genotype[i],evol.traitGenes[t]));
        mean /= metaPopSize();
        return mean;
    }

    double genotype(int p, int t)
    {
        double mean = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        		mean += Auxils.arrayMean(Auxils.arrayElements(genotype[i],evol.traitGenes[t]));
        mean /= popSize(p);
        return mean;
    }

    double genotypeVar(int t)
    {
        double mean = genotype(t);
        double var = 0;
        for (int i = 0; i < totSites; i++)
        	if (alive[i])
        		var += Math.pow(mean - Auxils.arrayMean(Auxils.arrayElements(genotype[i],evol.traitGenes[t])), 2);
        var /= metaPopSize();
        return var;
    }

    double genotypeVar(int p, int t)
    {
        double mean = genotype(p, t);
        double var = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        		var += Math.pow(mean - Auxils.arrayMean(Auxils.arrayElements(genotype[i],evol.traitGenes[t])), 2);
        var /= popSize(p);
        return var;
    }

    double phenotype(int t)
    {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
        	if (alive[i])
        		mean += traitFenotype[i][t];
        mean /= metaPopSize();
        return mean;
    }

    double phenotype(int p, int t)
    {
        double mean = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        		mean += traitFenotype[i][t];
        mean /= popSize(p);
        return mean;
    }

    double phenotypeVar(int t)
    {
        double mean = phenotype(t);
        double var = 0;
        for (int i = 0; i < totSites; i++)
        	if (alive[i])
        		var += Math.pow(mean - traitFenotype[i][t], 2);
        var /= metaPopSize();
        return var;
    }

    double phenotypeVar(int p, int t)
    {
        double mean = phenotype(p, t);
        double var = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        		var += Math.pow(mean - traitFenotype[i][t], 2);
        var /= popSize(p);
        return var;
    }

    double traitFitnessMax(int p, int t)
    {
    	double max = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i] && max < traitFitness[i][t])
        		max = traitFitness[i][t];
        return max;
    }

    double traitFitness(int p, int t)
    {
        double mean = 0;
        double max = traitFitnessMax(p, t);
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        		mean += traitFitness[i][t]/max;
        mean /= popSize(p);
        return mean;
    }

    double traitFitness(int p)
    {
    	double mean = 0;
    	for (int t = 0; t < comm.traits; t++)
    		mean += traitFitness(p, t);
    	mean /= comm.traits;
    	return mean;
   }

    double traitFitnessVar(int p)
    {
        double mean = traitFitness(p);
        double var = 0;
    	for (int t = 0; t < comm.traits; t++)
        		var += Math.pow(mean - traitFitness(p, t), 2);
        var /= popSize(p);
        return var;
    }

    double absFitness()
    {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
        	if (alive[i])
        		mean += fitness[i];
        mean /= metaPopSize();
        return mean;
    }

    double absFitness(int p)
    {
        double mean = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        		mean += fitness[i];
        mean /= popSize(p);
        return mean;
    }

 
    double fitness()
    {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
        	if (alive[i])
        		mean += fitness[i]/maxFitness[patch[i]];
        mean /= metaPopSize();
        return mean;
    }

    double fitness(int p)
    {
        double mean = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        		mean += fitness[i]/maxFitness[p];
//        		mean += fitness[i];
        mean /= popSize(p);
        return mean;
    }

    double fitnessVar(int p)
    {
        double mean = fitness(p);
        double var = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        		var += Math.pow(mean - fitness[i]/maxFitness[p], 2);
//        		var += Math.pow(mean - fitness[i], 2);
        var /= popSize(p);
        return var;
    }

	double fitnessMin(int p)
    {
    	double min = 1;
    	for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
    		if (alive[i])
    			if(fitness[i]/maxFitness[p] < min)
    				min = fitness[i]/maxFitness[p];
//    			if(fitness[i] < min)
//					min = fitness[i];
    	return min;
    }

    double fitnessMax(int p)
    {
    	double max = 0;
    	for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
    		if (alive[i])
    			if(fitness[i]/maxFitness[p] > max)
    				max = fitness[i]/maxFitness[p];
//    			if(fitness[i] > max)
//					max = fitness[i];
    	return max;
    }

    double load(int p)
    {
        double mean = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        		mean += 1 - fitness[i]/maxFitness[p];
//        		mean += 1 - fitness[i];
        mean /= popSize(p);
        return mean;
    }

    double loadVar(int p)
    {
        double mean = load(p);
        double var = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        		var += Math.pow(mean - (1 - fitness[i]/maxFitness[p]), 2);
//        		var += Math.pow(mean - (1 - fitness[i]), 2);
        var /= popSize(p);
        return var;
    }

    double selectionDiff(int p, int t)
    {
        double mean = 0;
        double sum = 0;
        double fitRel;
        double fitMean = fitness(p);
        double fenotpSd = Math.sqrt(phenotypeVar(p, t));
        double SDiff;
        
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        	{
        		fitRel = (fitness[i]/maxFitness[p])/fitMean;
        		mean += traitFenotype[i][t]*fitRel;
        		sum += fitRel;
        	}
        mean /= sum;
        if (fenotpSd == 0)
        	SDiff = 0;
        else
        	SDiff = Math.abs((mean - phenotype(p, t)))/fenotpSd;
        return SDiff;
    }

    double selectionDiff(int p)
    {
    	double mean = 0;
    	for (int t = 0; t < comm.traits; t++)
    		mean += selectionDiff(p, t);
    	mean /= comm.traits;
        return mean;
   }

    double selectionDiffVar(int p)
    {
        double mean = selectionDiff(p);
        double var = 0;
    	for (int t = 0; t < comm.traits; t++)
        		var += Math.pow(mean - selectionDiff(p, t), 2);
        var /= popSize(p);
        return var;
    }

    double pSex()
    {
        double mean = 0;
        for (int i = 0; i < totSites; i++)
        	if (alive[i])
        		mean += pSex[i];
        mean /= metaPopSize();
        return mean;
    }

    double pSex(int p)
    {
        double mean = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        		mean += pSex[i];
        mean /= popSize(p);
        return mean;
    }

    double pSexVar(int p)
    {
        double mean = pSex(p);
        double var = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i])
        		var += Math.pow(mean - pSex[i], 2);
        var /= popSize(p);
        return var;
    }
}


/* Ecological parameters/variables */
class Comm 
{
	String envType = "GLOBAL"; 
	String envSync = "NO"; 
	int envDims = 1;
    
    double[] patchDev = {1};
    double[] dimDev = {1};
    
    int traits = 2;
    double minEnv = 0.2;
    double maxEnv = 0.8;
    double sigmaE = 0.0;
    int microsites = 600;
    double d = 0.1;
    double[] demogrCost = {0.5};
    
    int gridSize = 2;
    int nbrPatches = gridSize*gridSize; 
    double pChange = 0.1;
    double[] envStep = {0.01};
    double[] pExt = {0};
    double[] dispRate = {0.01};
    double rho = 1;
    double pSex = 0;

    int [] traitDim;

    double[][] neighbours = new double[nbrPatches][nbrPatches];
    double[][] dispNeighbours = new double[nbrPatches][nbrPatches];

    void init()
    {
        nbrPatches = gridSize*gridSize;
        neighbours = new double[nbrPatches][nbrPatches];
        dispNeighbours = new double[nbrPatches][nbrPatches];
        calcDistNeighbours();
        
        traitDim = new int [traits];
        int dim = 0;
        for (int tr = 0; tr < traits; tr++) {
        	traitDim[tr] = dim++;
   			if (dim == envDims)
   				dim = 0;
        }
    }
    
    void calcDistNeighbours() 
    {
        for (int i = 0; i < gridSize; i++)
            for (int j = 0; j < gridSize; j++)
                for (int i2 = 0; i2 < gridSize; i2++)
                    for (int j2 = 0; j2 < gridSize; j2++) {
                       double dist = Math.sqrt(Math.pow(Math.min(Math.abs(i - i2),gridSize-Math.abs(i - i2)),2) + Math.pow(Math.min(Math.abs(j - j2),gridSize-Math.abs(j - j2)),2));
                       neighbours[j*gridSize+i][j2*gridSize+i2] = dist;
                    }
    }
    
    void calcDispNeighbours(int dr)
    {
        for (int i = 0; i < nbrPatches; i++) {
            for (int j = 0; j < nbrPatches; j++)
                dispNeighbours[i][j] = (i == j) ? 0 : (rho*Math.exp(-rho*neighbours[i][j]));
            double iSum = Auxils.arraySum(dispNeighbours[i]);
            for (int j = 0; j < nbrPatches; j++)
                dispNeighbours[i][j] = (i == j) ? (1 - dispRate[dr]) : (dispRate[dr]*dispNeighbours[i][j]/iSum);
        }
    }
}


/* Evolution parameters/variables */
class Evol 
{
    double omegaE = 0.02; 
    double divF = 1; 
    int traitLoci = 20;
    int lociPerTrait = traitLoci;
    int sexLoci = 10;
    int allLoci = traitLoci + sexLoci;
    double mutationRate = 1e-4;
    double mutationSize = 1;
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
    
    void init(Comm comm) 
    {
        divF = 2*Math.pow(Math.sqrt(comm.traits)*omegaE, 2);

        allLoci = traitLoci + sexLoci;
        
        lociPerTrait = traitLoci/comm.traits;

        allMother  = new int [allLoci];
        allFather  = new int [allLoci];
        allGenes = new int[2*allLoci];
        somMother  = new int [traitLoci];
        somFather  = new int [traitLoci];
        somGenes = new int[2*traitLoci];
        traitMother  = new int [comm.traits][lociPerTrait];
        traitFather  = new int [comm.traits][lociPerTrait];
        traitGenes = new int[comm.traits][2*lociPerTrait];
        sexMother  = new int [sexLoci];
        sexFather  = new int [sexLoci];
        sexGenes = new int[2*sexLoci];
      
/* somatic genes */
        for (int tr = 0; tr < comm.traits; tr++) {
        	for (int l = 0; l < lociPerTrait; l++) {
        		longPos = l + (tr * lociPerTrait);
        		traitMother[tr][l] =  longPos;
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
class Run 
{
    int runs = 1;
    int timeSteps = 10000;
    int printSteps = 100;
    int saveSteps = 1000;
    String fileName = "output_evolvingMetacommunity.csv";
}


/* initialize simulation run */
class Init 
{
	double pSex;

	double[][] environment;
    int[] N;
    double[][] genotype;
    
    public Init(Comm comm, Evol evol, int dd, int pd, int es) 
    {
		double dEnv;
    	environment = new double [comm.nbrPatches][comm.envDims];
        N = new int [comm.nbrPatches];
		genotype = new double [comm.nbrPatches][comm.traits];

		java.util.Arrays.fill(N, comm.microsites);

		if (comm.envType.equals("GLOBAL"))
		{
			for (int d = 0; d < comm.envDims; d++)
			{
				dEnv = comm.minEnv + (Auxils.random.nextDouble() * (comm.maxEnv - comm.minEnv));
				for (int p = 0; p < comm.nbrPatches; p++) 
					environment[p][d] = dEnv;
			}
		}
		else
		{
			for (int p = 0; p < comm.nbrPatches; p++)
				for (int d = 0; d < comm.envDims; d++) 
					environment[p][d] = comm.minEnv + (Auxils.random.nextDouble() * (comm.maxEnv - comm.minEnv));
		}

   		for (int p = 0; p < comm.nbrPatches; p++)
   		{
   			for (int tr = 0; tr < comm.traits; tr++) {
    			genotype[p][tr] = environment[p][comm.traitDim[tr]];
   			}
   		}
		pSex = comm.pSex >= 0 ? comm.pSex : Auxils.random.nextDouble();
    }


    
}


/* reading in parameter values from input file */
class Reader
{
	static void readInput(String fileName, Comm comm, Evol evol, Run run) throws IOException
	{
		try (BufferedReader input = new BufferedReader(new FileReader(fileName))) {
			String line;
			String[] words;
			int size;
			while ((line = input.readLine()) != null) {
				words = line.trim().split("\\s+");
				switch (words[0]) {
				case "ENVDIMS": comm.envDims = Integer.parseInt(words[1]);
				break;
				case "TRAITS": comm.traits = Integer.parseInt(words[1]);
				break;
				case "MINENV": comm.minEnv = Double.parseDouble(words[1]);
				break;
				case "MAXENV": comm.maxEnv = Double.parseDouble(words[1]);
				break;
				case "SIGMAE": comm.sigmaE = Double.parseDouble(words[1]);
				break;
				case "MICROSITES": comm.microsites = Integer.parseInt(words[1]);
				break;
				case "D": comm.d = Double.parseDouble(words[1]);
				break;
				case "PSEX": comm.pSex = Double.parseDouble(words[1]);
				break;
				case "COST": 
	            	size = Integer.parseInt(words[1]);
	            	comm.demogrCost = new double[size];
		            for (int i = 0; i < size; i++)
		                comm.demogrCost[i] = Double.parseDouble(words[2+i]);
				break;
				case "GRIDSIZE": comm.gridSize = Integer.parseInt(words[1]);
				break;
				case "ENVTYPE": comm.envType = words[1];
				break;
				case "ENVSYNC": comm.envSync = words[1];
				break;
				case "PCHANGE": comm.pChange = Double.parseDouble(words[1]);
				break;
				case "ENVSTEP": 
	            	size = Integer.parseInt(words[1]);
	            	comm.envStep = new double[size];
		            for (int i = 0; i < size; i++)
		                comm.envStep[i] = Double.parseDouble(words[2+i]);
				break;
				case "PEXT": 
	            	size = Integer.parseInt(words[1]);
	            	comm.pExt = new double[size];
		            for (int i = 0; i < size; i++)
		                comm.pExt[i] = Double.parseDouble(words[2+i]);
				break;
				case "M": 
	            	size = Integer.parseInt(words[1]);
	            	comm.dispRate = new double[size];
		            for (int i = 0; i < size; i++)
		                comm.dispRate[i] = Double.parseDouble(words[2+i]);
				break;
				case "RHO": comm.rho = Double.parseDouble(words[1]);
				break;

				case "OMEGAE": evol.omegaE = Double.parseDouble(words[1]);
				break;
				case "TRAITLOCI": evol.traitLoci = Integer.parseInt(words[1]);
				break;
				case "SEXLOCI": evol.sexLoci = Integer.parseInt(words[1]);
				break;
				case "MU": evol.mutationRate = Double.parseDouble(words[1]);
				break;
				case "MUSIZE": evol.mutationSize = Double.parseDouble(words[1]);
				break;
				case "SIGMAZ": evol.sigmaZ = Double.parseDouble(words[1]);
				break;

				case "RUNS": run.runs = Integer.parseInt(words[1]);
				break;
				case "TIMESTEPS": run.timeSteps = Integer.parseInt(words[1]);
				break;
				case "PRINTSTEPS": run.printSteps = Integer.parseInt(words[1]);
				break;
				case "SAVESTEPS": run.saveSteps = Integer.parseInt(words[1]);
				break;
				case "OUTPUT": run.fileName = words[1];
				break;
				}
			}

		}
	}
}


/* Auxiliary functions for array calculations */
class Auxils 
{
    static Random random = new Random();

    static void arrayShuffle(int[] array)
    {
        int index, temp;
        for (int i = array.length - 1; i > 0; i--)
        {
            index = random.nextInt(i + 1);
            temp = array[index];
            array[index] = array[i];
            array[i] = temp;
        }
    }

    static void arrayShuffle(double[] array)
    {
        int index;
        double temp;
        for (int i = array.length - 1; i > 0; i--)
        {
            index = random.nextInt(i + 1);
            temp = array[index];
            array[index] = array[i];
            array[i] = temp;
        }
    }
    
    static int[] arraySample(int n, int[] array)
    {
        int[] tempArr = array.clone();
        arrayShuffle(tempArr);
        return java.util.Arrays.copyOf(tempArr,n);
    }

    static double[] arraySample(int n, double[] array)
    {
        double[] tempArr = array.clone();
        arrayShuffle(tempArr);
        return java.util.Arrays.copyOf(tempArr,n);
    }

    //sampling with or without replacement
    static int[] arraySampleProb(int n, int[] array, double[] probs, boolean repl)
    {
    	int pos;
        double rand;
        int[] newElements;
        int[] newArr = new int[n];
        double[] cumProbs = probs.clone();
        arrayCumSum(cumProbs);
		arrayDiv(cumProbs, cumProbs[cumProbs.length-1]);
      		for (int i = 0; i < n; i++) {
            rand = random.nextDouble();
            pos = arraySearch(cumProbs, rand);
            newArr[i] = array[pos];
            if (!repl) {
            	newElements = arrayConcat(enumArray(0, pos-1), enumArray(pos+1, array.length-1));
            	array = arrayElements(array, newElements);
            	cumProbs = arrayElements(cumProbs, newElements);
                arrayCumSum(cumProbs);
                arrayDiv(cumProbs, cumProbs[cumProbs.length-1]);
            }
        }
        return newArr;
    }

    //sampling with or without replacement
    static double[] arraySampleProb(int n, double[] array, double[] probs, boolean repl)
    {
        int pos;
        double rand;
        int[] newElements;
        double[] newArr = new double[n];
        double[] cumProbs = probs.clone();
        arrayCumSum(cumProbs);
        arrayDiv(cumProbs, cumProbs[cumProbs.length-1]);
        for (int i = 0; i < n; i++) {
            rand = random.nextDouble();
            pos = arraySearch(cumProbs, rand);
            newArr[i] = array[pos];
            if (!repl) {
            	newElements = arrayConcat(enumArray(0, pos-1), enumArray(pos+1, array.length-1));
            	array = arrayElements(array, newElements);
            	cumProbs = arrayElements(cumProbs, newElements);
                arrayCumSum(cumProbs);
                arrayDiv(cumProbs, cumProbs[cumProbs.length-1]);
            }
        }
        return newArr;
    }

    static int randIntProb(int end, double[] probs)
    {
        int val;
        double rand;
        double[] cumProbs = java.util.Arrays.copyOf(probs, end);
        Auxils.arrayCumSum(cumProbs);
        Auxils.arrayDiv(cumProbs, cumProbs[cumProbs.length-1]);
        rand = random.nextDouble();
        val = arraySearch(cumProbs, rand);
        return val;
    }

    static int countDistinct(double[] arr, int n)
    { 
        // First sort the array so that all 
        // occurrences become consecutive 
    	java.util.Arrays.sort(arr); 
  
        // Traverse the sorted array 
        int res = 0; 
        for (int i = 0; i < n; i++)  
        { 
  
            // Move the index ahead while 
            // there are duplicates 
            while (i < n - 1 &&  
                    arr[i] == arr[i + 1]) 
            { 
                i++; 
            } 
            res++; 
        } 
        return res; 
    } 
    
    static int arraySearch(double[] array, double key)
    {
        int lo = 0;
        int hi = array.length - 1;
        int mid;
        if (key <= array[lo])
            return lo;
        else {
            while (lo < hi) {
                mid = lo + (hi - lo)/2;
                if (key <= array[mid])
                    hi = mid - 1;
                else if (key > array[mid])
                    lo = mid + 1;
            }
            if (key <= array[lo])
                return lo;
            else
                return ++lo;
        }
    }
    
    static int[] enumArray(int from, int to)
    {
        int[] newArr = new int[to-from+1];
        for (int i = 0; i < newArr.length; i++)
            newArr[i] = from++;
        return newArr;
    }
    
    static int[] arrayElements(int[] array, int[] pos)
    {
        int[] newArr = new int[pos.length];
        for (int i = 0; i < newArr.length; i++)
            newArr[i] = array[pos[i]];
        return newArr;
    }
    
    static boolean[] arrayElements(boolean[] array, int[] pos)
    {
        boolean[] newArr = new boolean[pos.length];
        for (int i = 0; i < newArr.length; i++)
            newArr[i] = array[pos[i]];
        return newArr;
    }
    
    static double[] arrayElements(double[] array, int[] pos)
    {
        double[] newArr = new double[pos.length];
        for (int i = 0; i < newArr.length; i++)
            newArr[i] = array[pos[i]];
        return newArr;
    }
    
    static double arrayMean(int[] array)
    {
        double mean = 0;
		for (int value : array) mean += value;
        mean /= array.length;
        return mean;
    }

    static double arrayMean(boolean[] array)
    {
        double mean = 0;
		for (boolean b : array)
			if (b)
				mean++;
        mean /= array.length;
        return mean;
    }

    static double arrayMean(double[] array)
    {
        double mean = 0;
		for (double v : array) mean += v;
        mean /= array.length;
        return mean;
    }

    static double arrayMean(int[] array, int end)
    {
        double mean = 0;
        for (int i = 0; i < end; i++)
            mean += array[i];
        mean /= end;
        return mean;
    }

    static double arrayMean(boolean[] array, int end)
    {
        double mean = 0;
        for (int i = 0; i < end; i++)
            if (array[i])
                mean++;
        mean /= end;
        return mean;
    }

    static double arrayMean(double[] array, int end)
    {
        double mean = 0;
        for (int i = 0; i < end; i++)
            mean += array[i];
        mean /= end;
        return mean;
    }

    static int arrayMax(int[] array)
    {
        int max = array[0];
        for (int i = 1; i < array.length; i++)
        	if (array[i] > max)
            max = array[i];
        return max;
    }

    static double arrayMax(double[] array)
    {
        double max = array[0];
        for (int i = 1; i < array.length; i++)
        	if (array[i] > max)
            max = array[i];
        return max;
    }

    static int arrayMin(int[] array)
    {
        int min = array[0];
        for (int i = 1; i < array.length; i++)
        	if (array[i] < min)
            min = array[i];
        return min;
    }

    static double arrayMin(double[] array)
    {
        double min = array[0];
        for (int i = 1; i < array.length; i++)
        	if (array[i] < min)
            min = array[i];
        return min;
    }

    static int arraySum(int[] array) 
    {
        int sum = 0;
		for (int value : array) sum += value;
        return sum;
    }

    static int arraySum(boolean[] array) 
    {
        int sum = 0;
		for (boolean b : array)
			if (b)
				sum++;
        return sum;
    }

    static double arraySum(double[] array) 
    {
        double sum = 0;
		for (double v : array) sum += v;
        return sum;
    }

    static int arraySum(int[] array, int end) 
    {
        int sum = 0;
        for (int i = 0; i < end; i++)
            sum += array[i];
        return sum;
    }

    static int arraySum(boolean[] array, int end) 
    {
        int sum = 0;
        for (int i = 0; i < end; i++)
            if (array[i])
                sum++;
        return sum;
    }

    static double arraySum(double[] array, int end) 
    {
        double sum = 0;
        for (int i = 0; i < end; i++)
            sum += array[i];
        return sum;
    }

    static void arrayCumSum(int[] array) 
    {
            for (int i = 1; i < array.length; i++)
                array[i] += array[i-1];
    }

    static void arrayCumSum(double[] array) 
    {
            for (int i = 1; i < array.length; i++)
                array[i] += array[i-1];
    }

    static void arrayAdd(int[] array, int a) 
    {
        for (int i = 0; i < array.length; i++)
            array[i] += a;
    }

    static void arrayAdd(double[] array, double a) 
    {
        for (int i = 0; i < array.length; i++)
            array[i] += a;
    }

    static void arrayMult(int[] array, int a) 
    {
        for (int i = 0; i < array.length; i++)
            array[i] *= a;
    }

    static void arrayMult(double[] array, double a) 
    {
        for (int i = 0; i < array.length; i++)
            array[i] *= a;
    }

    static void arrayDiv(int[] array, int a) 
    {
        for (int i = 0; i < array.length; i++)
            array[i] /= a;
    }

    static void arrayDiv(double[] array, double a) 
    {
        for (int i = 0; i < array.length; i++)
            array[i] /= a;
    }

    static int[] arrayConcat(int[] first, int[] second) 
    {
        int[] result = java.util.Arrays.copyOf(first, first.length + second.length);
        System.arraycopy(second, 0, result, first.length, second.length);
        return result;
    }    
    
    static double[] arrayConcat(double[] first, double[] second) 
    {
        double[] result = java.util.Arrays.copyOf(first, first.length + second.length);
        System.arraycopy(second, 0, result, first.length, second.length);
        return result;
    }    

    static double adjustToRange(double val, double min, double max) 
    {
    	double range = max - min;
    	int quot = (int) Math.floor((val - max)/range);
    	double rem = mod((val - max), range);
    	int minAdd = mod(quot, 2);
    	int maxAdd = 1 - minAdd;
		return minAdd*(min + rem) + maxAdd*(max - rem);
    }
    
    static int mod(int x, int y)
    {
        int result = x % y;
        return result < 0 ? result + y : result;
    }

    static double mod(double x, double y)
    {
        double result = x % y;
        return result < 0 ? result + y : result;
    }

    static double mod(int x, double y)
    {
        double result = x % y;
        return result < 0 ? result + y : result;
    }

    static double mod(double x, int y)
    {
        double result = x % y;
        return result < 0 ? result + y : result;
    }
}


