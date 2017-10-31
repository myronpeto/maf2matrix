package org.panda.resource.tcga;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

/**
 * Reads and serves MutSig results.
 *
 * @author Ozgun Babur
 */
public class MutSigReader
{
	public static Map<String, Double> readPValues(String dir)
	{
		return readValues(dir, true);
	}

	public static Map<String, Double> readQValues(String dir)
	{
		return readValues(dir, false);
	}

	public static boolean hasMutsig(String dir)
	{
		return new File(dir + File.separator + "scores-mutsig.txt").exists();
	}

	public static Map<String, Double> readValues(String dir, boolean p) { try
	{
		Map<String, Double> map = new HashMap<>();
		Scanner sc = new Scanner(new File(dir + File.separator + "scores-mutsig.txt"));
		sc.nextLine();
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			String gene = token[1];
			double pval = Double.parseDouble(token[token.length - (p ? 2 : 1)]);
			map.put(gene, pval);
		}
		return map;
	}
	catch (IOException e) {throw new RuntimeException(e);}}
}
