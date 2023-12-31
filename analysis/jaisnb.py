import numpy as np
import matplotlib.colors as col

jaisnb_rgb = np.array([[0.49019608, 0.,         0.49803922],
        [0.49803922, 0.,         0.52156863],
        [0.50588235, 0.,         0.54509804],
        [0.51372549, 0.,         0.56470588],
        [0.51764706, 0.,         0.58823529],
        [0.51764706, 0.,         0.60784314],
        [0.52156863, 0.,         0.62745098],
        [0.51764706, 0.,         0.64705882],
        [0.51764706, 0.,         0.6627451 ],
        [0.50980392, 0.,         0.68235294],
        [0.50588235, 0.,         0.69803922],
        [0.49803922, 0.,         0.71372549],
        [0.49019608, 0.,         0.72941176],
        [0.47843137, 0.,         0.74509804],
        [0.46666667, 0.,         0.76078431],
        [0.45490196, 0.,         0.77254902],
        [0.43921569, 0.,         0.78823529],
        [0.42352941, 0.,         0.8       ],
        [0.40784314, 0.,         0.81176471],
        [0.38823529, 0.,         0.82352941],
        [0.36862745, 0.,         0.83529412],
        [0.34901961, 0.,         0.84705882],
        [0.32941176, 0.,         0.85882353],
        [0.30980392, 0.,         0.86666667],
        [0.28627451, 0.,         0.87843137],
        [0.2627451 , 0.,         0.88627451],
        [0.23921569, 0.,         0.89411765],
        [0.21568627, 0.,         0.90196078],
        [0.19215686, 0.,         0.90980392],
        [0.16862745, 0.,         0.91764706],
        [0.14509804, 0.,         0.9254902 ],
        [0.11764706, 0.,         0.92941176],
        [0.09411765, 0.,         0.9372549 ],
        [0.06666667, 0.,         0.94117647],
        [0.03921569, 0.,         0.94901961],
        [0.01568627, 0.,         0.95294118],
        [0.,         0.00784314, 0.95686275],
        [0.,         0.03137255, 0.96078431],
        [0.,         0.05882353, 0.96470588],
        [0.,         0.08235294, 0.96862745],
        [0.,         0.10980392, 0.97254902],
        [0.,         0.13333333, 0.97647059],
        [0.,         0.16078431, 0.97647059],
        [0.,         0.18431373, 0.98039216],
        [0.,         0.20784314, 0.98431373],
        [0.,         0.23137255, 0.98431373],
        [0.,         0.25490196, 0.98823529],
        [0.,         0.27843137, 0.98823529],
        [0.,         0.30196078, 0.99215686],
        [0.,         0.3254902,  0.99215686],
        [0.,         0.34509804, 0.99215686],
        [0.,         0.36862745, 0.99215686],
        [0.,         0.38823529, 0.99607843],
        [0.,         0.40784314, 0.99607843],
        [0.,         0.42745098, 0.99607843],
        [0.,         0.44705882, 0.99607843],
        [0.,         0.46666667, 0.99607843],
        [0.,         0.48627451, 0.99607843],
        [0.,         0.50588235, 0.99607843],
        [0.,         0.52156863, 0.99607843],
        [0.,         0.54117647, 0.99607843],
        [0.00392157, 0.55686275, 0.99607843],
        [0.00392157, 0.57254902, 0.99607843],
        [0.00784314, 0.58823529, 0.99607843],
        [0.00784314, 0.60392157, 1.        ],
        [0.01176471, 0.61960784, 1.        ],
        [0.01568627, 0.63529412, 1.        ],
        [0.01960784, 0.65098039, 1.        ],
        [0.02352941, 0.66666667, 1.        ],
        [0.02745098, 0.67843137, 1.        ],
        [0.03137255, 0.69411765, 1.        ],
        [0.03529412, 0.70588235, 1.        ],
        [0.03921569, 0.71764706, 1.        ],
        [0.04705882, 0.73333333, 1.        ],
        [0.05098039, 0.74509804, 1.        ],
        [0.05882353, 0.75686275, 1.        ],
        [0.0627451,  0.76862745, 1.        ],
        [0.07058824, 0.78039216, 1.        ],
        [0.07843137, 0.79215686, 1.        ],
        [0.08235294, 0.8,        1.        ],
        [0.09019608, 0.81176471, 1.        ],
        [0.09803922, 0.82352941, 1.        ],
        [0.10980392, 0.83137255, 1.        ],
        [0.11764706, 0.84313725, 1.        ],
        [0.1254902,  0.85098039, 1.        ],
        [0.13333333, 0.85882353, 1.        ],
        [0.14509804, 0.87058824, 1.        ],
        [0.15294118, 0.87843137, 1.        ],
        [0.16470588, 0.88627451, 1.        ],
        [0.17647059, 0.89411765, 1.        ],
        [0.18431373, 0.90588235, 1.        ],
        [0.19607843, 0.91372549, 1.        ],
        [0.20784314, 0.92156863, 1.        ],
        [0.21960784, 0.92941176, 1.        ],
        [0.23137255, 0.9372549,  1.        ],
        [0.24705882, 0.94509804, 1.        ],
        [0.25882353, 0.95294118, 1.        ],
        [0.27058824, 0.96078431, 1.        ],
        [0.28627451, 0.96862745, 1.        ],
        [0.29803922, 0.97254902, 1.        ],
        [0.31372549, 0.98039216, 1.        ],
        [0.32941176, 0.98823529, 1.        ],
        [0.34509804, 0.99607843, 1.        ],
        [0.35686275, 1.,         0.99607843],
        [0.37647059, 1.,         0.98823529],
        [0.39215686, 1.,         0.98431373],
        [0.40784314, 1.,         0.97647059],
        [0.42352941, 1.,         0.97254902],
        [0.44313725, 1.,         0.96862745],
        [0.45882353, 1.,         0.96078431],
        [0.47843137, 1.,         0.95686275],
        [0.49803922, 1.,         0.95294118],
        [0.51764706, 1.,         0.94901961],
        [0.5372549,  1.,         0.94901961],
        [0.55686275, 1.,         0.94509804],
        [0.58039216, 1.,         0.94117647],
        [0.60392157, 1.,         0.94117647],
        [0.62352941, 1.,         0.9372549 ],
        [0.64705882, 1.,         0.9372549 ],
        [0.6745098,  1.,         0.9372549 ],
        [0.69803922, 1.,         0.94117647],
        [0.7254902,  1.,         0.94117647],
        [0.75294118, 1.,         0.94509804],
        [0.78431373, 1.,         0.94901961],
        [0.81960784, 1.,         0.95294118],
        [0.85490196, 1.,         0.96078431],
        [0.89803922, 1.,         0.97254902],
        [0.95294118, 1.,         0.98431373],
        [0.95686275, 1.,         0.95294118],
        [1.,1.,1.],
        [1.,1.,1.],
        [0.90980392, 1.,         0.89803922],
        [0.8745098,  1.,         0.85490196],
        [0.84705882, 1.,         0.81960784],
        [0.82352941, 1.,         0.78431373],
        [0.80392157, 1.,         0.75294118],
        [0.78431373, 1.,         0.7254902 ],
        [0.76862745, 1.,         0.69803922],
        [0.75294118, 1.,         0.6745098 ],
        [0.7372549,  1.,         0.64705882],
        [0.72941176, 1.,         0.62352941],
        [0.71764706, 1.,         0.60392157],
        [0.70980392, 1.,         0.58039216],
        [0.70196078, 1.,         0.55686275],
        [0.69803922, 1.,         0.5372549 ],
        [0.69019608, 1.,         0.51764706],
        [0.68627451, 1.,         0.49803922],
        [0.68627451, 1.,         0.47843137],
        [0.68235294, 1.,         0.45882353],
        [0.68235294, 1.,         0.44313725],
        [0.68627451, 1.,         0.42352941],
        [0.68627451, 1.,         0.40784314],
        [0.68627451, 1.,         0.39215686],
        [0.69019608, 1.,         0.37647059],
        [0.69411765, 1.,         0.35686275],
        [0.69803922, 1.,         0.34509804],
        [0.70588235, 1.,         0.32941176],
        [0.70980392, 1.,         0.31372549],
        [0.71764706, 1.,         0.29803922],
        [0.7254902,  1.,         0.28627451],
        [0.73333333, 1.,         0.27058824],
        [0.74117647, 1.,         0.25882353],
        [0.74901961, 1.,         0.24705882],
        [0.75686275, 1.,         0.23137255],
        [0.76862745, 1.,         0.21960784],
        [0.77647059, 1.,         0.20784314],
        [0.78823529, 1.,         0.19607843],
        [0.8,        1.,         0.18431373],
        [0.80784314, 1.,         0.17647059],
        [0.81960784, 1.,         0.16470588],
        [0.83137255, 1.,         0.15294118],
        [0.84313725, 1.,         0.14509804],
        [0.85490196, 1.,         0.13333333],
        [0.86666667, 1.,         0.1254902 ],
        [0.87843137, 1.,         0.11764706],
        [0.89019608, 1.,         0.10980392],
        [0.90196078, 1.,         0.09803922],
        [0.91372549, 1.,         0.09019608],
        [0.9254902,  1.,         0.08235294],
        [0.9372549,  1.,         0.07843137],
        [0.94901961, 1.,         0.07058824],
        [0.96078431, 1.,         0.0627451 ],
        [0.97254902, 1.,         0.05882353],
        [0.98431373, 1.,         0.05098039],
        [0.99607843, 1.,         0.04705882],
        [1.,         0.98823529, 0.03921569],
        [1.,         0.97647059, 0.03529412],
        [1.,         0.96470588, 0.03137255],
        [1.,         0.95294118, 0.02745098],
        [1.,         0.94117647, 0.02352941],
        [1.,         0.92941176, 0.01960784],
        [1.,         0.91764706, 0.01568627],
        [1.,         0.90588235, 0.01176471],
        [1.,         0.89411765, 0.00784314],
        [0.99607843, 0.88235294, 0.00784314],
        [0.99607843, 0.8745098,  0.00392157],
        [0.99607843, 0.8627451,  0.00392157],
        [0.99607843, 0.85098039, 0.        ],
        [0.99607843, 0.83921569, 0.        ],
        [0.99607843, 0.82745098, 0.        ],
        [0.99607843, 0.81568627, 0.        ],
        [0.99607843, 0.8,        0.        ],
        [0.99607843, 0.78823529, 0.        ],
        [0.99607843, 0.77647059, 0.        ],
        [0.99607843, 0.76470588, 0.        ],
        [0.99607843, 0.75294118, 0.        ],
        [0.99215686, 0.7372549,  0.        ],
        [0.99215686, 0.7254902,  0.        ],
        [0.99215686, 0.70980392, 0.        ],
        [0.99215686, 0.69803922, 0.        ],
        [0.98823529, 0.68235294, 0.        ],
        [0.98823529, 0.66666667, 0.        ],
        [0.98431373, 0.65098039, 0.        ],
        [0.98431373, 0.63529412, 0.        ],
        [0.98039216, 0.61960784, 0.        ],
        [0.97647059, 0.60392157, 0.        ],
        [0.97647059, 0.58823529, 0.        ],
        [0.97254902, 0.56862745, 0.        ],
        [0.96862745, 0.55294118, 0.        ],
        [0.96470588, 0.53333333, 0.        ],
        [0.96078431, 0.51372549, 0.        ],
        [0.95686275, 0.49411765, 0.        ],
        [0.95294118, 0.4745098,  0.        ],
        [0.94901961, 0.45490196, 0.        ],
        [0.94117647, 0.43529412, 0.        ],
        [0.9372549,  0.41568627, 0.        ],
        [0.92941176, 0.39215686, 0.        ],
        [0.9254902,  0.37254902, 0.        ],
        [0.91764706, 0.35294118, 0.        ],
        [0.90980392, 0.32941176, 0.        ],
        [0.90196078, 0.30588235, 0.        ],
        [0.89411765, 0.28627451, 0.        ],
        [0.88627451, 0.2627451,  0.        ],
        [0.87843137, 0.23921569, 0.        ],
        [0.86666667, 0.21568627, 0.        ],
        [0.85882353, 0.19215686, 0.        ],
        [0.84705882, 0.17254902, 0.        ],
        [0.83529412, 0.14901961, 0.        ],
        [0.82352941, 0.1254902,  0.        ],
        [0.81176471, 0.10196078, 0.        ],
        [0.8,        0.08235294, 0.        ],
        [0.78823529, 0.05882353, 0.        ],
        [0.77254902, 0.03921569, 0.        ],
        [0.76078431, 0.01568627, 0.        ],
        [0.74509804, 0.,         0.        ],
        [0.72941176, 0.,         0.01960784],
        [0.71372549, 0.,         0.03921569],
        [0.69803922, 0.,         0.05490196],
        [0.68235294, 0.,         0.0745098 ],
        [0.6627451,  0.,         0.09019608],
        [0.64705882, 0.,         0.10588235],
        [0.62745098, 0.,         0.12156863],
        [0.60784314, 0.,         0.13333333],
        [0.58823529, 0.,         0.14509804],
        [0.56470588, 0.,         0.15686275],
        [0.54509804, 0.,         0.16862745],
        [0.52156863, 0.,         0.17647059],
        [0.49803922, 0.,         0.18039216]])

jaisnb = col.ListedColormap(jaisnb_rgb, name='jaisnb')

