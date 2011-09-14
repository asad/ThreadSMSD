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
package thread;

import helper.SMSDThread;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class SMSDTests {

    public static void main(String[] args) {
        try {
            String[] smiles = new String[]{
                "O=C(OCCN(C)C)C=C",
                "O=C(O)C3N1C(=O)CC1S(=O)(=O)C3(C)(CN2N=NC=C2)",
                "O=NN2C=C(C=1C=CC=CC=12)CCO",
                "O=NN(C=1C=CC=CC=1)C",
                "O=C=NC1=CC=CC(N=C=O)=C1C",
                "O=C1C=C(OC2=CC(OCC(=O)OCC)=CC=C12)C3=CC=CC=C3",
                "N1=C(NC(=C1C=2C=CC=CC=2)C=3C=CC=CC=3)C4=CC=C(OC)C=C4",
                "O=C5C=1C=CC=CC=1N(C3=C5(C(O)=CC=2OC(CC=23)C4(OC4)(C)))C",
                "C=1C(=C(C=C(C=1Cl)Cl)Cl)Cl",
                "O=C(NC1=CC=CC(=C1)Br)C=2C=C([N+](=O)[O-])SC=2",
                "C1=CC=CC=C1", 
                "C1=CC=NC=C1",
                "CCO",
                "C1=CC=NC=C1",
                "C1=CC=C(C=C1)C1=CC=CC=C1",
                "C1=CC=C(C=C1)[N+]1=CC=CC=C1",
                "C1=CC2=C3C(C=CC=N3=CC=C2)=C1",
                "C1=CC2=CC=C(C=C2C=C1)C1=C2C=CC=CC2=CC=C1",
                "C1=CC2=C(C=C1)C=CC=C2"
            };

            ExecutorService executorService = //Executors.newCachedThreadPool();
                    Executors.newFixedThreadPool(smiles.length);

            IAtomContainer[] compounds = new IAtomContainer[smiles.length];
            SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
            for (int i = 0; i < smiles.length; i++) {
                compounds[i] = sp.parseSmiles(smiles[i]);
            }


            List<Future<List<String>>> futureList = new ArrayList<Future<List<String>>>();
            for (int i = 0; i < compounds.length; i++) {
                if (compounds[i] != null) {
                    for (int j = i + 1; j < compounds.length; j++) {
                        if (compounds[i] != null) {
                            try {

                                Callable<List<String>> getMCSS = new SMSDThread(
                                        compounds[i],
                                        compounds[j]);
                                Future<List<String>> future = executorService.submit(getMCSS);
                                futureList.add(future);
                            } catch (Exception e) {
                                e.printStackTrace();
                            }
                        }
                    }
                }
            }
            //Display the results....
            for (Iterator<Future<List<String>>> it = futureList.iterator(); it.hasNext();) {
                //getting the future result (all smiles)
                Future<List<String>> future = it.next();

                List<String> smileList = future.get();
                for (String smileString : smileList) {
                    System.out.println("Smile of substructure is : " + smileString);
                }
            }
            System.out.println("\nDone\n");
            executorService.shutdown();
        } catch (Exception e) {
            e.printStackTrace();
        }

        System.out.println("\nFINALLY Done\n");
    }
}
