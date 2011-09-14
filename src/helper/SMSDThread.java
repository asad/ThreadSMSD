/*
 * Copyright (C) 2011 Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package helper;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.Callable;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.interfaces.Algorithm;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class SMSDThread implements Callable<List<String>> {

    private final IAtomContainer compound1;
    private final IAtomContainer compound2;

    /**
     * Constructor. Initializes with the URI of the compounds
     * @param compound1 query
     * @param compound2 target
     */
    public SMSDThread(
            IAtomContainer compound1,
            IAtomContainer compound2) {
        //now you need to remove the hydrogen...
        this.compound1 = AtomContainerManipulator.removeHydrogens(compound1);
        this.compound2 = AtomContainerManipulator.removeHydrogens(compound2);
    }

    /**
     * The call function implemented from the <code>Callable</code> interface
     * 
     * @return string of the smile retrieved
     */
    @Override
    public synchronized List<String> call() {
        try {
            List<IAtomContainer> mcssList = Collections.synchronizedList(findMCSSList());
            List<IAtomContainer> finalIAtomContainerList = Collections.synchronizedList(new ArrayList<IAtomContainer>());
            List<String> fragmentSMILEsList = Collections.synchronizedList(new ArrayList<String>());
            Iterator<IAtomContainer> subgraphIterator = mcssList.iterator();

            SmilesGenerator sg = new SmilesGenerator();
            while (subgraphIterator.hasNext()) {
                IAtomContainer overlap = subgraphIterator.next();
                String smile = sg.createSMILES(new Molecule(overlap));
                if (overlap.getAtomCount() > 0) {
                    fragmentSMILEsList.add(smile);
                    finalIAtomContainerList.add(overlap);
                }
                //System.out.println("Smile found : <" + smile + "> atom count=" + overlap.getAtomCount());
            }

            return fragmentSMILEsList;
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }

    }

    public synchronized List<IAtomContainer> findMCSSList() throws Exception {
        List<IAtomContainer> mcssList = Collections.synchronizedList(new ArrayList<IAtomContainer>());
        try {

            boolean bondSensitive = true;
            boolean ringMatcher = true;
            boolean stereoMatch = true;
            boolean fragmentMinimization = true;
            boolean energyMinimization = true;

            Isomorphism comparison = new Isomorphism(getQueryCompound(), getTargetCompound(), Algorithm.VFLibMCS, bondSensitive, ringMatcher);
            comparison.setChemFilters(stereoMatch, fragmentMinimization, energyMinimization);

            int count_final_sol = 1;
            //System.out.println("Output of the final Mappings: ");

            if (comparison.getAllAtomMapping() != null) {
                for (AtomAtomMapping final_solution : comparison.getAllAtomMapping()) {
                    Collection<IAtom> queryMatches = final_solution.getMappings().keySet();
                    IAtomContainer fragment = getMatchedSubgraph(getQueryCompound(), queryMatches);
                    mcssList.add(fragment);

                    System.out.println("Stereo Match: " + comparison.getStereoScore(count_final_sol - 1));
                    System.out.println("Stereo different: " + comparison.isStereoMisMatch());
                    System.out.println("Fragment Size: " + comparison.getFragmentSize(count_final_sol - 1));
                    System.out.println("Tanimoto Similarity Score: " + comparison.getTanimotoSimilarity());
                    System.out.println("Tanimoto Euclidean Distance: " + comparison.getEuclideanDistance());
                    count_final_sol++;
                }
            }
            return mcssList;
        } catch (Exception ex) {
            ex.printStackTrace();
            return mcssList;
        }
    }

    public synchronized List<IAtomContainer> findSubstructureList() throws Exception {
        List<IAtomContainer> mcssList = Collections.synchronizedList(new ArrayList<IAtomContainer>());
        try {

            boolean bondSensitive = true;
            boolean ringMatcher = true;
            boolean stereoMatch = true;
            boolean fragmentMinimization = true;
            boolean energyMinimization = true;

            Substructure comparison = new Substructure(getQueryCompound(), getTargetCompound(), bondSensitive, ringMatcher, true);
            comparison.setChemFilters(stereoMatch, fragmentMinimization, energyMinimization);

            int count_final_sol = 1;
            //System.out.println("Output of the final Mappings: ");

            if (comparison.getAllAtomMapping() != null) {
                for (AtomAtomMapping final_solution : comparison.getAllAtomMapping()) {
                    Collection<IAtom> queryMatches = final_solution.getMappings().keySet();
                    IAtomContainer fragment = getMatchedSubgraph(getQueryCompound(), queryMatches);
                    mcssList.add(fragment);

                    System.out.println("Stereo Match: " + comparison.getStereoScore(count_final_sol - 1));
                    System.out.println("Stereo different: " + comparison.isStereoMisMatch());
                    System.out.println("Fragment Size: " + comparison.getFragmentSize(count_final_sol - 1));
                    System.out.println("Tanimoto Similarity Score: " + comparison.getTanimotoSimilarity());
                    System.out.println("Tanimoto Euclidean Distance: " + comparison.getEuclideanDistance());
                    count_final_sol++;
                }
            }
            return mcssList;
        } catch (Exception ex) {
            ex.printStackTrace();
            return mcssList;
        }
    }

    /**
     * This function will returns a subgraph extracted from the
     * source or target molecules
     * @param container source/target container
     * @param matches   source/target mapping
     * @return mapped subgraph/substructure
     */
    public synchronized IAtomContainer getMatchedSubgraph(IAtomContainer container, Collection<IAtom> matches) {
        IAtomContainer needle = container.getBuilder().newInstance(IAtomContainer.class, container);
        List<IAtom> atomListToBeRemoved = new ArrayList<IAtom>();
        for (IAtom containerAtom : container.atoms()) {
            boolean discardAtom = true;
            for (IAtom matchedAtom : matches) {
                if (containerAtom == matchedAtom) {
                    discardAtom = false;
                }
            }
            if (discardAtom) {
                int index = container.getAtomNumber(containerAtom);
                atomListToBeRemoved.add(needle.getAtom(index));
            }
        }
        for (IAtom removeAtom : atomListToBeRemoved) {
            needle.removeAtomAndConnectedElectronContainers(removeAtom);
        }
        atomListToBeRemoved.clear();
        return needle;
    }

    /**
     * @return the Query Compound
     */
    public synchronized IAtomContainer getQueryCompound() {
        return compound1;
    }

    /**
     * @return the Target Compound
     */
    public synchronized IAtomContainer getTargetCompound() {
        return compound2;
    }
}
