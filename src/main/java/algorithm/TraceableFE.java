package algorithm;

import it.unisa.dia.gas.jpbc.*;
import it.unisa.dia.gas.plaf.jpbc.pairing.PairingFactory;
import it.unisa.dia.gas.plaf.jpbc.pairing.a.TypeACurveGenerator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.*;

public class TraceableFE {

    public static final String paramFilePath = "params.properties";
    public static final String propertyFilePath = "properties.properties";

    class Params {
        Pairing pairing;
        Field G;
        Element g1;
        Element g2;
        BigInteger p;
        Element egg;
        Gamma gamma;

        public Params(Pairing pairing, Field G, Element g1, Element g2, BigInteger p, Element egg, Gamma gamma) {
            this.pairing = pairing;
            this.G = G;
            this.g1 = g1;
            this.g2 = g2;
            this.p = p;
            this.egg = egg;
            this.gamma = gamma;
        }
    }

    class PK {
        Params params;
        List<Element> hs;
        List<Element> bs;

        public PK(Params params, List<Element> hs, List<Element> bs) {
            this.params = params;
            this.hs = hs;
            this.bs = bs;
        }
    }

    class MSK {

        List<Element> ss;
        List<Element> ts;

        public MSK(List<Element> ss, List<Element> ts) {
            this.ss = ss;
            this.ts = ts;
        }
    }

    class SKX {
        public SKX(Element tk) {
            this.tk = tk;
        }

        Element tk;
    }

    class Gamma {

        List<Theta> gamma;

        public Gamma(List<Theta> gamma) {
            this.gamma = gamma;
        }
    }

    class Theta {
        List<Element> theta;

        public Theta(List<Element> theta) {
            this.theta = theta;
        }
    }

    class SetupParams {
        public SetupParams(PK PK, MSK msk) {
            this.PK = PK;
            this.msk = msk;
        }

        PK PK;
        MSK msk;
    }

    class CT0 {
        Element ct0;

        public CT0(Element ct0) {
            this.ct0 = ct0;
        }
    }

    class CTI {
        List<Element> cti;

        public CTI(List<Element> cti) {
            this.cti = cti;
        }
    }

    class CT {
        List<Element> hg;
        List<Element> br;

        public CT(List<Element> hg, List<Element> br) {
            this.hg = hg;
            this.br = br;
        }
    }

    class DecryptResult {
        public DecryptResult(Element decryptResult) {
            this.decryptResult = decryptResult;
        }

        Element decryptResult;
    }

    public void curveInit(int lambda) throws IOException {
        // curve 256bit, limit field 512bit
        // A E A1 F four curves
        PairingParametersGenerator generator = new TypeACurveGenerator(lambda, lambda * 2);
        PairingParameters parameters = generator.generate();
        BufferedWriter writer = new BufferedWriter(new FileWriter(paramFilePath));
        writer.write(parameters.toString());
        writer.close();
    }

    public SetupParams setup(int lambda, int k, int n) throws IOException {
        // init curve params
        File paramFile = new File(paramFilePath);
        if (!paramFile.exists() || paramFile.length() == 0) {
            curveInit(lambda);
        }

        // get params
        Pairing pairing = PairingFactory.getPairing(paramFilePath);
        Field G1 = pairing.getG1();
        Field G2 = pairing.getG2();
        BigInteger q = G1.getOrder();
        Element g1 = G1.newRandomElement().getImmutable();
        Element g2 = G2.newRandomElement().getImmutable();

        // init t and b
        List<Element> ts = new ArrayList<>();
        List<Element> bs = new ArrayList<>();
        for (int i = 0; i < k; i++) {
            ts.add(pairing.getZr().newRandomElement().getImmutable());
            bs.add(g1.powZn(ts.get(i)).getImmutable());
        }

        // init G h and s
        Element egg = pairing.pairing(g1, g2).getImmutable();
        List<Element> hs = new ArrayList<>();
        List<Element> ss = new ArrayList<>();
        for (int i = 0; i < k; i++) {
            ss.add(pairing.getZr().newRandomElement().getImmutable());
            hs.add(egg.powZn(ss.get(i)).getImmutable());
        }

        // init gamma
        List<TraceableFE.Theta> thetas = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            List<Element> theta = new ArrayList<>();
            for (int j = 0; j < k; j++) {
                theta.add(pairing.getZr().newRandomElement().getImmutable());
            }
            thetas.add(new Theta(theta));
        }
        TraceableFE.Gamma gamma = new Gamma(thetas);

        // init params
        Params params = new Params(pairing, G1, g1, g2, q, egg, gamma);

        // init mpk and msk
        PK pk = new PK(params, hs, bs);
        MSK msk = new MSK(ss, ts);

        // init setup params
        SetupParams setupParams = new SetupParams(pk, msk);

        return setupParams;
    }

    public SKX extract(Params params, int id, MSK msk, List<Element> x) {

        // get thetaID
        Theta thetaID = params.gamma.gamma.get(id);

        // get s and t
        List<Element> ss = msk.ss;
        List<Element> ts = msk.ts;

        // compute intermedia result
        Element sx = params.pairing.getZr().newZeroElement();
        Element tTheta = params.pairing.getZr().newZeroElement();
        for (int i = 0; i < x.size(); i++) {
            sx = sx.add(ss.get(i).mul(x.get(i)));
            tTheta = tTheta.add(ts.get(i).mul(thetaID.theta.get(i)));
        }

        // compute tk
        Element tk = sx.div(tTheta).getImmutable();
        return new SKX(tk);
    }

    public CT encrypt(PK pk, List<Element> y) {
        Params params = pk.params;

        // random choose
        Element r = params.pairing.getZr().newRandomElement().getImmutable();

        // compute ct
        List<Element> hg = new ArrayList<>();
        List<Element> br = new ArrayList<>();
        for (int i = 0; i < y.size(); i++) {
            Element hr = pk.hs.get(i).powZn(r);
            Element gy = params.egg.powZn(y.get(i));
            hg.add(hr.mul(gy).getImmutable());
            br.add(pk.bs.get(i).powZn(r));
        }
        return new CT(hg, br);
    }

    public DecryptResult decrypt(Params params, CT ct, SKX SKX, List<Element> x, int id) {
        // get thetaid
        Theta thetaID = params.gamma.gamma.get(id);

        // compute upper part and lower part
        Element decryptResult = params.pairing.getGT().newOneElement();
        Element brTheta = params.pairing.getG1().newOneElement();
        for (int i = 0; i < x.size(); i++) {
            decryptResult = decryptResult.mul(ct.hg.get(i).powZn(x.get(i)));
            brTheta = brTheta.mul(ct.br.get(i).powZn(thetaID.theta.get(i)));
        }
        Element g2tk = params.g2.powZn(SKX.tk).getImmutable();

        // pair get the lower part
        Element ebrThetag2tk = params.pairing.pairing(brTheta, g2tk);

        // div to get result
        decryptResult = decryptResult.div(ebrThetag2tk);
        return new DecryptResult(decryptResult);
    }
}