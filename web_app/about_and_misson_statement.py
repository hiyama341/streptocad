from dash import dcc, html
from styling import text_style

about_streptocad_and_mission_statement = html.Div([
    dcc.Markdown(
        """
        ### **StreptoCAD's Mission**

        
        - Our mission is to accelerate biological research and innovation, enabling you to design and build Streptomyces strains faster and with greater precision through automated workflows.


        ### **Scaling Biology for Faster Learning**

        
        - As new challenges emerge, we leverage synthetic biology to engineer solutions. From therapeutic microbes that combat diseases to biofactories producing valuable compounds, and ecosystems designed to mitigate environmental impacts, our goal is to push the boundaries of what's possible in biological engineering with automation and software.


        ### **Join the Community**

        
        - Become a part of the StreptoCAD community and collaborate with fellow innovators. Share your workflows, gain insights from others, and contribute to a collective effort in advancing synthetic biology. Create and customize your own workflows to fit your unique research needs, and help shape the future of biological engineering.


        ### **Happy bioengineering!**
        """,
        style={**text_style, "marginBottom": "20px", "lineHeight": "2", "marginTop": "30px"}
    ),
    html.Img(
        src="assets/intro_fig.png",  
        style={"width": "50%", "height": "auto", "display": "block", "marginLeft": "0", "marginRight": "auto"}  
    )
])
